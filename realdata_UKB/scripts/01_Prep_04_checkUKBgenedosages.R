#' ---
#' title: "Check UKB genetic data (candidate SNPs)"
#' subtitle: "Longitudinal MVMR"
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
#' Here, I want to load the genetic data. 
#' 
#' In my PLINK2 call, I already filter for imputation quality (mach r2 >0.8) and MAF (>0.01). In the SNP selection script, I already selected the variants I want to test and for the samples with longitudinal data.
#' 
#' Here, I just load the pgen data and store in an format for R for the later GAMLSS runs. Then, I check per chromosome the pairwise LD.
#' 
#' # Init ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile.R")
.libPaths()

#' # Load genetic data ####
#' ***
UKB_genotypes_filtered = gsub("~","../../../",UKB_genotypes_filtered)
pvar1 = NewPvar(paste0(UKB_genotypes_filtered, 'UKB_filtered.pvar'))
pgen = NewPgen(paste0(UKB_genotypes_filtered,'UKB_filtered.pgen'), pvar=pvar1)
pvar = fread(paste0(UKB_genotypes_filtered,'UKB_filtered.pvar'))
psam = fread(paste0(UKB_genotypes_filtered,'UKB_filtered.psam'))

load("../results/01_Prep_03_SNPList.RData")
table(pvar$ID == SNPList$rsID)
setnames(pvar,"#CHROM","CHR")

myNRs = 1:dim(pvar)[1]
geno_mat <- ReadList(pgen, variant_subset = myNRs , meanimpute=F)
dim(geno_mat)
colnames(geno_mat) = pvar[,ID]
rownames(geno_mat) = psam$IID

# filt samples
load(paste0(UKB_phenotypes_filtered,"/01_Prep_01_BL_FU_GP_merged_filtered.RData"))
filt = psam$IID %in% myTab_long$ID
table(filt)
head(psam)
setnames(psam, "#FID","FID")

psam = psam[filt,]
geno_mat = geno_mat[filt,]

pvar[,EAF := colSums(geno_mat)/(2*dim(psam)[1])]
pvar[,MAF := EAF]
pvar[EAF>0.5,MAF := 1-EAF]
pvar[,table(MAF<0.01)]
pvar_AF = fread(paste0(UKB_genotypes_filtered, '/UKB_merged_AF.afreq'))
pvar_AF = pvar_AF[ID %in% pvar$ID]
stopifnot(pvar$ID == pvar_AF$ID)
pvar[,EAF_fullUKB := pvar_AF$ALT_FREQS]
plot(pvar$EAF,pvar$EAF_fullUKB)
pvar[,EAF_fullUKB := NULL]

save(pvar, psam, geno_mat, file = paste0(UKB_genotypes_filtered,"UKB_merged.RData"))

#' # LD checks ####
#' ***
pvar[,comment_LD := "LD OK"]
test = pvar[,.N,CHR]
test = test[N>1,]
myCHRs = test$CHR

LDTab = foreach(i = 1:length(myCHRs))%do%{
  #i=1
  myCHR = myCHRs[i]
  filt = pvar$CHR == myCHR
  dumMatrix = geno_mat[,filt]
  CorTab = cor(dumMatrix)^2
  heatmap(CorTab,Rowv = NA,Colv = NA, main =paste0("Chromosome ",myCHRs[i]))
  CorTab2 = as.data.table(CorTab)
  CorTab2[,SNP1 := rownames(CorTab)]
  CorTab_long = melt(CorTab2,id.vars="SNP1",measure.vars=rownames(CorTab),variable.name="SNP2")
  CorTab_long = CorTab_long[SNP1 != SNP2,]
  CorTab2_long = CorTab_long[value>0.05,]
  pvar[ID %in% CorTab2_long$SNP1 | ID %in% CorTab2_long$SNP2, comment_LD := "LD >0.05 with at least 1 SNP"]
  CorTab_long
}
LDTab = rbindlist(LDTab)
dev.off()
table(pvar$comment_LD)

#' Okay, there are many SNPs with no LD proxy, and some SNP which have at least one LD proxy with LD r2>0.05. 
#'
#' I will decide which SNP to use when I have the SNP associations. Then I can pick the best of these correlated variants. 
#' 
#' # Save data ####
#' ***
save(LDTab, file = paste0("../results/01_Prep_04_LD.RData"))

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
