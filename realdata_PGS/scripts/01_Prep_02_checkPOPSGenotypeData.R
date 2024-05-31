#' ---
#' title: "Check POPS genetic data (PGS SNPs)"
#' subtitle: "Longitudinal MVMR in POPS"
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
#' Here I want to load the provided SNP data and make some initial checks (without any phenotype information)
#' 
#' - check AF between POPS and reference
#' - check HWE (imputed SNP usually will not be HWE controlled as imputation takes care of that, but I just want to check anyway)
#' - check LD (can I use all SNPs or just some SNPs?)
#'
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_HPC.R")
source("../../helperfunctions/HWETest.R")
.libPaths()

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("2024-","24-",tag)
tag = gsub("-","",tag)

#' # Load data ####
#' ***
#' 
#' - PLINK2 files: pgen (genotype matrix), pvar (SNP info), psam (sample info)
#' - phenotype files: taken from script 02
#' - SNP list: taken from script 01
#' 
#' ## PLINK2 data ####
pvar1 = NewPvar(paste0(POPS_SNP_data_PGS, '.pvar'))
pgen = NewPgen(paste0(POPS_SNP_data_PGS, '.pgen'), pvar=pvar1)
pvar = fread(paste0(POPS_SNP_data_PGS, '.pvar'))
psam = fread(paste0(POPS_SNP_data_PGS, '.psam'))

geno_mat <- ReadList(pgen, variant_subset = c(1:dim(pvar)[1]) , meanimpute=F)
dim(geno_mat)

#' some modifcations...
head(psam)
setnames(psam, "#FID","FID")

head(pvar)
setnames(pvar, "#CHROM","CHR")
dummy = pvar$INFO
dummy2 = unlist(strsplit(dummy,";"))
filt_AF = grepl("AF=",dummy2) & !grepl("MAF=",dummy2)
filt_MAF = grepl("MAF=",dummy2)
filt_R2 = grepl("R2=",dummy2) & !grepl("ER2=",dummy2)
filt_ER2 = grepl("ER2=",dummy2)
filt_type = grepl("IMPUTED",dummy2) | grepl("TYPED",dummy2)

dummy_MAF = dummy2[filt_MAF]
dummy_MAF = gsub("MAF=","",dummy_MAF)
pvar[, MAF:= as.numeric(dummy_MAF)]

dummy_AF = dummy2[filt_AF]
dummy_AF = gsub("AF=","",dummy_AF)
pvar[, AF:= as.numeric(dummy_AF)]

dummy_R2 = dummy2[filt_R2]
dummy_R2 = gsub("R2=","",dummy_R2)
pvar[, R2:= as.numeric(dummy_R2)]

dummy_type = dummy2[filt_type]
pvar[, type:= dummy_type]

dummy_ER2 = dummy2[filt_ER2]
dummy_ER2 = gsub("ER2=","",dummy_ER2)
pvar[type == "TYPED", ER2:= as.numeric(dummy_ER2)]

class(geno_mat)
colnames(geno_mat) = pvar$ID
rownames(geno_mat) = psam$FID

#' ## SNP list ####
myFiles = list.files("../results/",pattern = "01_Prep_01")
myFiles = myFiles[grep("RData",myFiles)]
myFiles
loaded = load(paste0("../results/",myFiles))
SNPList = copy(get(loaded))

table(is.element(SNPList$SNP,pvar$ID))

SNPList = SNPList[SNP %in% pvar$ID, ]
SNPList = SNPList[, ID := SNP]
matched = match(pvar$ID,SNPList$ID)
SNPList = SNPList[matched,]

#' # Checks ####
#' ***
#' ## Allele frequency ####
stopifnot(pvar$ID == SNPList$ID)
plot(pvar$MAF,SNPList$MAF)
plot(pvar$AF,SNPList$ALT_Frq)

table(pvar$ALT == SNPList$`ALT(1)`)
table(pvar$REF == SNPList$`REF(0)`)
filt = pvar$ALT == SNPList$`ALT(1)`
plot(pvar$AF[filt],SNPList$ALT_Frq[filt])
plot(pvar$AF[!filt],SNPList$ALT_Frq[!filt])

#' Okay, there are a handful of SNPs which might have coding issues: their alleles should not be the same, but the allele frequencies are. This will be a problem when comparing the effect directions (planned sanity check). I will flag those SNPs for later. 
#' 
SNPList[,EAF2 := ALT_Frq]
SNPList[!filt, EAF2 := 1-EAF2]
AF_dif = abs(pvar$AF - SNPList$EAF2)
table(AF_dif>0.1)
filt = AF_dif >= 0.01
max(AF_dif[!filt])
min(AF_dif[filt])
plot(pvar$AF[!filt],SNPList$EAF2[!filt])
points(pvar$AF[filt],SNPList$EAF2[filt],col="red")
pvar[,comment_AF := "AF OK"]
pvar[filt,comment_AF := "AF switched?"]
SNPList[,comment_AF := "AF OK"]
SNPList[filt,comment_AF := "AF switched?"]

#' ## Imputation quality ####
#' How good is the imputation of the selected SNPs?
hist(pvar$R2)
table(pvar$R2<0.8)
pvar[R2<0.8,c(1:5,8:11)]

#' Okay, there are 318 SNPs with low imputation quality. For the moment I will keep them, I can filter for R2 later.  
#' 
pvar[,comment_R2 := "R2 OK"]
pvar[R2<0.8,comment_R2 := "R2 low"]
SNPList[,comment_R2 := pvar$comment_R2]

#' ## Hardy-Weinberg #### 
#' Are the SNPs violating HWE? (no filter as I use imputed data, but would be nice to know)
HWETab = foreach(i=1:dim(pvar)[1])%do%{
  #i=1
  mySNP = geno_mat[,i]
  filt_AA = mySNP <= 0.5
  filt_AB = mySNP > 0.5 & mySNP < 1.5
  filt_BB = mySNP >= 1.5
  
  mySNP[filt_AA] = 0
  mySNP[filt_AB] = 1
  mySNP[filt_BB] = 2
  
  table1 = table(mySNP)
  names(table1)[names(table1)==0] = "AA"
  names(table1)[names(table1)==1] = "AB"
  names(table1)[names(table1)==2] = "BB"
  
  AA = table1[names(table1)=="AA"]
  AB = table1[names(table1)=="AB"]
  BB = table1[names(table1)=="BB"]
  
  if(length(AA)==0) AA = 0
  if(length(AB)==0) AB = 0
  if(length(BB)==0) BB = 0
  
  test = HWETest(AA = AA, 
                 AB = AB,
                 BB = BB)
  test[,NR := i]
  test
}
HWETab = rbindlist(HWETab)
table(HWETab$Pvalue<=0.05)
table(HWETab$Pvalue<=0.05/dim(pvar)[1])
table(HWETab$Pvalue<=1e-6)

#' Okay, there are some SNPs violating HWE. I will highlight them as well. 
#' 
pvar[,comment_HWE := "HWE OK"]
pvar[HWETab$Pvalue<1e-6,comment_HWE := "HWE pval sig"]
SNPList[,comment_HWE := pvar$comment_HWE]

#' Create a filter for high quality SNPs
pvar[,highQC := T]
pvar[!grepl("OK",comment_AF),highQC := F]
pvar[!grepl("OK",comment_R2),highQC := F]
pvar[!grepl("OK",comment_HWE),highQC := F]
table(pvar$highQC)

#' Okay, there are 5517 SNPs with high quality (same AF as in GWAS catalog, high R2, no sig HWE violation). I will test all 5861, but might decide later to only use the 5517 highQC SNPs for instrument selection. 
#' 
#' ## Pairwise LD ####
#' are the SNPs truly not in LD? (per chromosome, testing simple correlation)
pvar[,comment_LD := "LD OK"]
myCHRs = unique(pvar$CHR)
LDTab = foreach(i = 1:length(myCHRs))%do%{
  #i=1
  myCHR = myCHRs[i]
  filt = pvar$CHR == myCHR
  dumMatrix = geno_mat[,filt]
  CorTab = cor(dumMatrix)^2
  heatmap(CorTab,Rowv = NA,Colv = NA, main =paste0("Chromosome ",i))
  CorTab2 = as.data.table(CorTab)
  CorTab2[,SNP1 := rownames(CorTab)]
  CorTab_long = melt(CorTab2,id.vars="SNP1",measure.vars=rownames(CorTab),variable.name="SNP2")
  CorTab_long = CorTab_long[SNP1 != SNP2,]
  CorTab2_long = CorTab_long[value>0.05,]
  pvar[ID %in% CorTab2_long$SNP1 | ID %in% CorTab2_long$SNP2, comment_LD := "LD >0.05 with at least 1 SNP"]
  CorTab_long
}
LDTab = rbindlist(LDTab)

table(pvar$comment_LD)

#' Okay, there are 3335 SNPs with no LD proxy, and 2526 SNP which have at least one LD proxy with LD r2>0.05. I will filter for independent SNPs after getting the SNP statistics, as I want to do priority pruning depending on association p-values. 
#' 
#' # Save data ####
#' ***
save(pvar, psam, geno_mat, file = paste0(POPS_phenotypes,"/01_Prep_02_SNPData_",tag,".RData"))

save(HWETab, file= paste0("../results/01_Prep_02_HWE_",tag,".RData"))
save(LDTab, file = paste0("../results/01_Prep_02_LD_",tag,".RData"))
save(SNPList, file = paste0("../results/01_Prep_02_SNPList_",tag,".RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
