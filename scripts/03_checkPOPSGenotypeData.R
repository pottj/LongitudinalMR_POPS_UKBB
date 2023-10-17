#' ---
#' title: "Check POPS genetic data"
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
#' Here I want to load the provided SNP data and make some initial checks
#' 
#' - are all IDs with phenotype data available?
#' - is the babys sex in the phenotype data equal to the genetic sex?
#' - what are the allele frequency in POPS compared to the GWAS catalog data? (compare coding alleles) 
#' - how good is the imputation of the selected SNPs?
#' - are the SNPs violating HWE? (no filter as I use imputated data, but would be nice to know)
#' - are the SNPs truly not in LD? (per chromosome, testing simple correlation)
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_HPC.R")
.libPaths()

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("2023-","23-",tag)
tag = gsub("-","",tag)

#' # Load data ####
#' ***
#' 
#' - PLINK2 files: pgen (genotype matrix), pvar (SNP info), psam (sample info)
#' - phenotype files: taken from script 02
#' - SNP list: taken from script 01
#' 
#' ## PLINK2 data ####
pvar1 = NewPvar(paste0(POPS_SNP_data, '.pvar'))
pgen = NewPgen(paste0(POPS_SNP_data, '.pgen'), pvar=pvar1)
pvar = fread(paste0(POPS_SNP_data, '.pvar'))
psam = fread(paste0(POPS_SNP_data, '.psam'))

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

#' ## Phenotype data ####
load("../data/IndividualLevelData/02_LongitudinalExposure.RData")
load("../data/IndividualLevelData/02_Outcome.RData")

#' ## SNP list ####
SNPList = fread("../results/01_SNPList_BirthWeight_231017.txt")

#' # Checks ####
#' ***
#' ## IDs ####
#' - are all IDs with phenotype data available?
table(is.element(myTab_X$POPSID,psam$FID))
myTab_X[!is.element(POPSID,psam$FID)]

table(is.element(myTab_Y$POPSID,psam$FID))
myTab_Y[!is.element(POPSID,psam$FID)]

table(is.element(psam$FID,myTab_Y$POPSID))

#' There are two IDs with no matching genotype data, and 478 IDS with no mathing phenotype data. I will filter for IDs with both phnotype and genotype data. In addition, I will match the data (easier to read and analyze later)
#' 
goodSamples = myTab_Y[is.element(POPSID,psam$FID),POPSID]
myTab_Y = myTab_Y[is.element(POPSID,psam$FID),]
myTab_X = myTab_X[is.element(POPSID,psam$FID),]
matched = match(myTab_Y$POPSID,psam$FID)
psam = psam[matched,]
table(psam$FID == myTab_Y$POPSID)
geno_mat = geno_mat[matched,]

#' ## Sex ####
#' - is the babys sex in the phenotype data equal to the genetic sex?
table(psam$SEX,myTab_Y$pn_sex)
filt = psam$SEX == 2 & myTab_Y$pn_sex == "MALE"
psam[filt,]
myTab_Y[filt,]

#' There is one sample with mismatching sex. I will exclude this one, but check with Jasmine if this is a known mismatch or if it a case of chromosomal aberation. I will ask how they typically handled it. 
#' 
myTab_Y = myTab_Y[!filt,]
myTab_X = myTab_X[POPSID %in% myTab_Y$POPSID,]
psam = psam[!filt,]
geno_mat = geno_mat[!filt,]

#' ## Allele frequency ####
#' - what are the allele frequency in POPS compared to the GWAS catalog data? (compare coding alleles) 
pvar[,dumID := paste(CHR, POS, sep="_")]
SNPList[,dumID := paste(chr,pos, sep="_")]
table(duplicated(pvar$dumID))
table(is.element(pvar$dumID,SNPList$dumID))
matched = match(pvar$dumID,SNPList$dumID)
SNPList = SNPList[matched,]
table(pvar$dumID == SNPList$dumID)
plot(pvar$AF,SNPList$EAF)
SNPList[,MAF := EAF]
SNPList[EAF>0.5,MAF := 1-EAF]
plot(pvar$MAF,SNPList$MAF)

table(pvar$ALT == SNPList$EA)
filt = pvar$ALT == SNPList$EA
plot(pvar$AF[filt],SNPList$EAF[filt])
plot(pvar$AF[!filt],SNPList$EAF[!filt])

#' Okay, there are a handful of SNPs which might have coding issues: their MAFs are similar in the SNP list and genotype data, but the EAFs are not. This will be a problem when comparing the effect directions (planned sanity check). I have more faith in the actual genotype data of POPS than the GWAS catalog data. It is highly likely that the effect allele was wrongly noted there. I will flag those SNPs for later. 
#' 
SNPList[,EAF2 := EAF]
SNPList[!filt, EAF2 := 1-EAF2]
AF_dif = abs(pvar$AF - SNPList$EAF2)
table(AF_dif>0.1)
filt = AF_dif >= 0.1
max(AF_dif[!filt])
min(AF_dif[filt])
plot(pvar$AF[!filt],SNPList$EAF2[!filt])
points(pvar$AF[filt],SNPList$EAF2[filt],col="red")
pvar[,comment_AF := "AF OK"]
pvar[filt,comment_AF := "AF switched?"]
SNPList[,comment_AF := "AF OK"]
SNPList[filt,comment_AF := "AF switched?"]

#' ## Imputation quality ####
#' - how good is the imputation of the selected SNPs?
hist(pvar$R2)
table(pvar$R2<0.8)
pvar[R2<0.8,c(1:5,8:11,14)]

#' Okay, there are 7 SNPs with low imputation quality. For the moment I will keep them, I can filter for R2 later.  
#' 
pvar[,comment_R2 := "R2 OK"]
pvar[R2<0.8,comment_R2 := "R2 low"]
SNPList[,comment_R2 := pvar$comment_R2]

#' ## Hardy-Weinberg #### 
#' - are the SNPs violating HWE? (no filter as I use imputated data, but would be nice to know)
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

#' Okay, there are some SNPs violating HWE. I will highlight them as well. 
#' 
pvar[,comment_HWE := "HWE OK"]
pvar[HWETab$Pvalue<0.05/dim(pvar)[1],comment_HWE := "HWE pval sig"]
SNPList[,comment_HWE := pvar$comment_HWE]

#' Create a filter for high quality SNPs
pvar[,highQC := T]
pvar[!grepl("OK",comment_AF),highQC := F]
pvar[!grepl("OK",comment_R2),highQC := F]
pvar[!grepl("OK",comment_HWE),highQC := F]
table(pvar$highQC)

#' Okay, there are 332 SNPs with high quality (same AF as in GWAS catalog, high R2, no sig HWE violation). I will test all 365, but might decide later to only use the 332 highQC SNPs for instrument selection. 
#' 
#' ## Pairwise LD ####
#' - are the SNPs truly not in LD? (per chromosome, testing simple correlation)
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

#' Okay, there are 169 SNPs with no LD proxy, and 196 SNP which have at least one LD proxy with LD r2>0.05. I will filter for independent SNPs after getting the SNP statistics, as I want to do priority pruning depending on association p-values. 
#' 
#' # Save data ####
#' ***
save(pvar, psam, geno_mat, file = paste0("../data/IndividualLevelData/03_POPS_SNPData_filtered_",tag,".RData"))
save(myTab_X,file = paste0("../data/IndividualLevelData/03_LongitudinalExposure_filtered_",tag,".RData"))
save(myTab_Y,file = paste0("../data/IndividualLevelData/03_Outcome_filtered_",tag,".RData"))

save(HWETab, file= paste0("../results/03_HWE_",tag,".RData"))
save(LDTab, file = paste0("../results/03_LD_",tag,".RData"))
save(SNPList, file = paste0("../results/03_SNPList_",tag,".RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
