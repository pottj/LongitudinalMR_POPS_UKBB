#' ---
#' title: "Check candidate SNPs from GWAS Catalog"
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
#' In the real data application of my longitudinal MVMR analyses, I want to use SNPs affecting heart rate (HR).
#' 
#' Here, I load the top SNPs for HR from the GWAS Catalog (downloaded 2024-06-27): 
#' 
#' - EFO_0004326: heart rate                                      --> 9184 (EFO changed)
#' - EFO_0004351: resting heart rate                              --> 4351
#' - EFO_0008003: heart rate variability measurement              --> 9258 (EFO changed)
#' - EFO_0009184: heart rate response to exercise                 --> 9184
#' - EFO_0009185: heart rate response to recovery post exercise   --> 9185
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_HPC.R")
.libPaths()

#' # Load data ####
#' ***
myFiles = list.files(path = pathData, pattern = "gwas-association-downloaded_2024-06-27")

dumTab = foreach(i = 1:length(myFiles))%do%{
  #i=1
  data = fread(paste0(pathData,myFiles[i]))
  data = data[PVALUE_MLOG >= -log10(5e-8),]
  data = data[!is.na(SNPS),]
  data = data[grepl("rs",SNPS),]
  data[,source := gsub("gwas-association-downloaded_2024-06-27-","",myFiles[i])]
  data[,source := gsub(".tsv","",source)]
  data[,source := gsub("-withChildTraits","",source)]
  
  data
}
myTab = rbindlist(dumTab)
myTab[,RAF := as.numeric(`RISK ALLELE FREQUENCY`)]
myTab[,MAF := RAF]
myTab[RAF>0.5,MAF := 1-RAF]

myTab[,table(duplicated(SNPS))]
myTab[,table(is.na(RAF),duplicated(SNPS))]

myTab3 = copy(myTab)
setorder(myTab3,-PVALUE_MLOG)
myTab3 = myTab3[!duplicated(SNPS),]
myTab3[,table(CHR_ID)]

#' # Save data #### 
#' ***
save(myTab3,file = paste0(UKB_phenotypes_filtered,"/01_Prep_03_SNPList_HR.RData"))
save(myTab3,file = "../results/01_Prep_03_SNPList_HR.RData")
write.table(myTab3$SNPS,file = paste0(UKB_phenotypes_filtered,"/01_Prep_03_SNPList_HR.txt"), 
            col.names = F, row.names = F, quote = F)
write.table(myTab3$SNPS,file = "../results/01_Prep_03_SNPList_HR.txt", 
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

myCHR = unique(myTab3$CHR_ID)
myCHR = myCHR[!is.na(myCHR)]
myCHR = myCHR[order(myCHR)]

dumTab = foreach(i = 1:length(myCHR))%do%{
  #i=1
  call2 = paste0("plink2", 
                   " --bgen ",UKB_SNP_data,"/ukb22828_c",myCHR[i],"_b0_v3.bgen",
                   " 'ref-last'", 
                   " --sample ",UKB_SNP_data,"ukb22828_c",myCHR[i],"_b0_v3_s487160.sample",
                   " --chr ", myCHR[i],
                   " --keep-fam ",UKB_phenotypes_filtered,"/01_Prep_02_SampleList_HR.txt",
                   " --extract ",UKB_phenotypes_filtered,"/01_Prep_03_SNPList_HR.txt", 
                   " --mach-r2-filter 0.8 2",
                   " --maf 0.01", 
                   " --threads 20",
                   " --export bgen-1.2 bits=8 id-delim='-'",
                   " --out ",UKB_genotypes_filtered,"/temp/UKB_HR_chr",myCHR[i])
  print(call2)
  out_filename = paste0(UKB_genotypes_filtered,"/temp/UKB_HR_chr",myCHR[i],".bgen")
  out_filename
}

call3 = c("cat-bgen -clobber -g")
               
for(i in 1:length(myCHR)){
  #i=1
  call3 = paste(call3,dumTab[[i]])
}

call3 = paste0(call3, " -og ",UKB_genotypes_filtered, "/temp/UKB_HR_merged.bgen")
print(call3)

call4 = paste0("plink2",
               " --bgen ",UKB_genotypes_filtered, "/temp/UKB_HR_merged.bgen",
               " 'ref-last'",
               " --sample ",UKB_genotypes_filtered, "/temp/UKB_HR_chr1.sample",
               " --make-pgen --out ",UKB_genotypes_filtered,"/UKB_HR_merged")
print(call4)

call5 = paste0("rm -rf ",UKB_genotypes_filtered,"/temp")
print(call5)

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

