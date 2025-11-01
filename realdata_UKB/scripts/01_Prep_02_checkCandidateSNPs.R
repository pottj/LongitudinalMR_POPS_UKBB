#' ---
#' title: "Check candidate SNPs from GLGC"
#' subtitle: "Longitudinal MR"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document:
#'     toc: true
#'     number_sections: true
#'     toc_float: true
#'     code_folding: show
#' ---
#'
#' # Introduction ####
#' ***
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile.R")
.libPaths()

#' # Load data ####
#' ***
myFiles = c(GLGC_TC_comb,GLGC_TC_females,GLGC_TC_males)
myFiles

dumTab1 = foreach(i = 1:length(myFiles))%do%{
  #i=1
  mySetting = gsub(".meta.singlevar.results.gz","",myFiles[i])
  mySetting = gsub(".*_","",mySetting)
  message("Working on ",mySetting)
  
  data = fread(myFiles[i])
  message("     Initial SNP number: ",dim(data)[1])
  
  # some filtering
  data = data[pvalue_neg_log10_GC >=-log10(5e-8),]
  data = data[!is.na(rsID)]
  data = data[POOLED_ALT_AF >= 0.01]
  data = data[POOLED_ALT_AF <= 0.99]
  data[,dumID := paste(CHROM,POS_b37,sep=":")]
  dupPos = data[duplicated(dumID)]
  data = data[!is.element(dumID,dupPos$dumID)]
  data[,dumID := NULL]
  message("     SNP number after filtering: ",dim(data)[1])
  
  # add column for sex-setting
  data[,setting := mySetting]
  
  # return data table
  data
}
TC_data = rbindlist(dumTab1)
dim(TC_data)
setorder(TC_data,CHROM,POS_b37,setting)
length(unique(TC_data$rsID))

#' # Save data #### 
#' ***
save(TC_data, file="../results/01_Prep_02_SNPList_GLGC.RData")
#load("../results/01_Prep_02_SNPList_GLGC.RData")
write.table(unique(TC_data$rsID),file = paste0(UKB_genotypes_filtered,"/01_Prep_02_SNPList_GLGC.txt"), 
            col.names = F, row.names = F, quote = F)

#' # Create PLINK2 calls ####
#' ***
#' I generate the calls here, and then copy-paste them into a slurm script. 
#' 
#' 1) create temp folder
#' 2) create bgen per chromosome (stored in temp)
#' 3) merge bgen files into one file (stored in temp)
#' 4) create pgen merged file (stored **not** in temp)
#' 5) remove temp folder (do not store all the chromosome-wise data)
#' 
call1 = paste0("mkdir ",UKB_genotypes_filtered,"/temp")
print(call1)

myCHR = unique(TC_data$CHROM)
myCHR = myCHR[!is.na(myCHR)]
myCHR = myCHR[order(myCHR)]

dumTab = foreach(i = 1:length(myCHR))%do%{
  #i=1
  call2 = paste0("plink2", 
                 " --bgen ",UKB_SNP_data,"/ukb22828_c",myCHR[i],"_b0_v3.bgen",
                 " 'ref-last'", 
                 " --sample ",UKB_SNP_data,"ukb22828_c",myCHR[i],"_b0_v3_s487160.sample",
                 " --chr ", myCHR[i],
                 " --keep-fam ",UKB_phenotypes_filtered,"/01_Prep_01_UKB_filtered_samples.txt",
                 " --extract ",UKB_genotypes_filtered,"/01_Prep_02_SNPList_GLGC.txt", 
                 " --mach-r2-filter 0.8 2",
                 " --maf 0.01", 
                 " --threads 20",
                 " --export bgen-1.2 bits=8 id-delim='-'",
                 " --out ",UKB_genotypes_filtered,"/temp/UKB_chr",myCHR[i])
  print(call2)
  out_filename = paste0(UKB_genotypes_filtered,"/temp/UKB_chr",myCHR[i],".bgen")
  out_filename
}

call3 = c("cat-bgen -clobber -g")

for(i in 1:length(myCHR)){
  #i=1
  call3 = paste(call3,dumTab[[i]])
}

call3 = paste0(call3, " -og ",UKB_genotypes_filtered, "/temp/UKB_merged.bgen")
print(call3)

call4 = paste0("plink2",
               " --bgen ",UKB_genotypes_filtered, "/temp/UKB_merged.bgen",
               " 'ref-last'",
               " --sample ",UKB_genotypes_filtered, "/temp/UKB_chr1.sample",
               " --make-pgen --out ",UKB_genotypes_filtered,"/UKB_merged")
print(call4)

call5 = paste0("rm -rf ",UKB_genotypes_filtered,"/temp")
print(call5)

call6 = paste0("plink2",
               " --pfile ",UKB_genotypes_filtered, "/UKB_merged",
               " --freq --out ",UKB_genotypes_filtered,"/UKB_merged_AF")
print(call6)

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
