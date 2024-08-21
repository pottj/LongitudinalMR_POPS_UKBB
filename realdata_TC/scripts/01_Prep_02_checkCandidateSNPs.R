#' ---
#' title: "Check candidate SNPs from TrajGWAS publication"
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

source("../../SourceFile_HPC.R")
.libPaths()

#' # Load data ####
#' ***
myFiles = list.files(path = paste0(pathData,"../TrajGWAS/"),pattern = "ukbbgen")
myFiles
myFiles = myFiles[grepl("totchol",myFiles)]

data = fread(paste0(pathData,"../TrajGWAS/",myFiles))
data[,table(betapval<5e-8,taupval<5e-8)]
data = data[betapval<5e-8 | taupval<5e-8,]
data[,table(maf<0.01,taupval<5e-8)]
data = data[maf>=0.01,]
data[,pheno := "TC"]
data

length(unique(data$snpid))
setorder(data,chr,pos)

# remove duplicates (triallelic SNPs)
data[,dumID := paste0(chr,":",pos)]
data[,table(duplicated(dumID))]
dups = data[duplicated(dumID),dumID]
data[dumID %in% dups]
data = data[!is.element(dumID, dups)]
data = data[grepl("rs",snpid),]
data[,table(chr)]

#' # Save data #### 
#' ***
save(data,file = paste0(UKB_phenotypes_filtered,"/01_Prep_02_SNPList_TC.RData"))
save(data, file="../results/01_Prep_02_SNPList_TC.RData")
write.table(data$snpid,file = paste0(UKB_phenotypes_filtered,"/01_Prep_02_SNPList_TC.txt"), 
            col.names = F, row.names = F, quote = F)
write.table(data$snpid,file = "../results/01_Prep_02_SNPList_TC.txt", 
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

myCHR = unique(data$chr)
myCHR = myCHR[!is.na(myCHR)]
myCHR = myCHR[order(myCHR)]

dumTab = foreach(i = 1:length(myCHR))%do%{
  #i=1
  call2 = paste0("plink2", 
                   " --bgen ",UKB_SNP_data,"/ukb22828_c",myCHR[i],"_b0_v3.bgen",
                   " 'ref-last'", 
                   " --sample ",UKB_SNP_data,"ukb22828_c",myCHR[i],"_b0_v3_s487160.sample",
                   " --chr ", myCHR[i],
                   " --keep-fam ",UKB_phenotypes_filtered,"/01_Prep_01_SampleList_TC.txt",
                   " --extract ",UKB_phenotypes_filtered,"/01_Prep_02_SNPList_TC.txt", 
                   " --mach-r2-filter 0.8 2",
                   " --maf 0.01", 
                   " --threads 20",
                   " --export bgen-1.2 bits=8 id-delim='-'",
                   " --out ",UKB_genotypes_filtered,"/temp/UKB_TC_chr",myCHR[i])
  print(call2)
  out_filename = paste0(UKB_genotypes_filtered,"/temp/UKB_TC_chr",myCHR[i],".bgen")
  out_filename
}

call3 = c("cat-bgen -clobber -g")
               
for(i in 1:length(myCHR)){
  #i=1
  call3 = paste(call3,dumTab[[i]])
}

call3 = paste0(call3, " -og ",UKB_genotypes_filtered, "/temp/UKB_TC_merged.bgen")
print(call3)

call4 = paste0("plink2",
               " --bgen ",UKB_genotypes_filtered, "/temp/UKB_TC_merged.bgen",
               " 'ref-last'",
               " --sample ",UKB_genotypes_filtered, "/temp/UKB_TC_chr1.sample",
               " --make-pgen --out ",UKB_genotypes_filtered,"/UKB_TC_merged")
print(call4)

call5 = paste0("rm -rf ",UKB_genotypes_filtered,"/temp")
print(call5)

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

