#' ---
#' title: "Get data for UKB CAD associations"
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
load("../temp/01_Prep_01_UKBB_filtered_Ancestry_Kinship_Meds.RData")
write.table(myTab$ID,file = paste0(UKB_phenotypes_filtered,"/01_Prep_09_SampleList_CAD.txt"),
            col.names = F, row.names = F, quote = F)

myPhenoFile1 = copy(myTab)
myPhenoFile1 = myPhenoFile1[,c(1,1,22,22,22,2)]
names(myPhenoFile1)[1:5] = c("FID","IID","CAD_combined","CAD_men","CAD_women")
myPhenoFile1[sex==1, CAD_women := NA]
myPhenoFile1[sex==0, CAD_men := NA]
myPhenoFile1[,sex:=NULL]
myPhenoFile1[,CAD_combined:= 1+ CAD_combined]
myPhenoFile1[,CAD_women:= 1+ CAD_women]
myPhenoFile1[,CAD_men:= 1+ CAD_men]
write.table(myPhenoFile1,file=paste0(UKB_phenotypes_filtered,"/01_Prep_09_Pheno_CAD.txt"),
            sep=" ",col.names = T,row.names = F,quote=F)

myCovarFile1 = copy(myTab)
myCovarFile1 = myCovarFile1[,c(1,1,2,6,20,7:16)]
names(myCovarFile1)[1:2] = c("FID","IID")
myCovarFile1[sex==0,sex:=2]
write.table(myCovarFile1,file=paste0(UKB_phenotypes_filtered,"/01_Prep_09_Covar_CAD.txt"),
            sep=" ",col.names = T,row.names = F,quote=F)

load("../results/01_Prep_02_SNPList_TC.RData")

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
                   " --keep-fam ",UKB_phenotypes_filtered,"/01_Prep_09_SampleList_CAD.txt",
                   " --extract ",UKB_phenotypes_filtered,"/01_Prep_02_SNPList_TC.txt", 
                   " --mach-r2-filter 0.8 2",
                   " --maf 0.01", 
                   " --threads 20",
                   " --export bgen-1.2 bits=8 id-delim='-'",
                   " --out ",UKB_genotypes_filtered,"/temp/UKB_CAD_chr",myCHR[i])
  print(call2)
  out_filename = paste0(UKB_genotypes_filtered,"/temp/UKB_CAD_chr",myCHR[i],".bgen")
  out_filename
}

call3 = c("cat-bgen -clobber -g")
               
for(i in 1:length(myCHR)){
  #i=1
  call3 = paste(call3,dumTab[[i]])
}

call3 = paste0(call3, " -og ",UKB_genotypes_filtered, "/temp/UKB_CAD_merged.bgen")
print(call3)

call4 = paste0("plink2",
               " --bgen ",UKB_genotypes_filtered, "/temp/UKB_CAD_merged.bgen",
               " 'ref-last'",
               " --sample ",UKB_genotypes_filtered, "/temp/UKB_CAD_chr1.sample",
               " --make-pgen --out ",UKB_genotypes_filtered,"/UKB_CAD_merged")
print(call4)

call5 = paste0("rm -rf ",UKB_genotypes_filtered,"/temp")
print(call5)

call6 = paste0("plink2", 
                 " --pfile ",UKB_genotypes_filtered,"/UKB_CAD_merged",
                 " --glm hide-covar firth-fallback",
                 " cols=chrom,pos,ref,alt,firth,test,nobs,machr2,a1freq,a1freqcc,a1countcc,orbeta,se,ci,tz,p",
                 " --pheno ",UKB_phenotypes_filtered,"/01_Prep_09_Pheno_CAD.txt", 
                 " --covar ",UKB_phenotypes_filtered,"/01_Prep_09_Covar_CAD.txt",
                 " --threads 20",
                 " --out ",UKB_genotypes_filtered,"/01_Prep_09_CAD")
call6

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

