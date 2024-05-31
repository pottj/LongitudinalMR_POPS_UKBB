#' ---
#' title: "Check POPS genetic data (EGG SNPs)"
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
POPS_SNP_data_EGG2 = gsub("~","../../..",POPS_SNP_data_EGG)
pvar1 = NewPvar(paste0(POPS_SNP_data_EGG2, '.pvar'))
pgen = NewPgen(paste0(POPS_SNP_data_EGG2, '.pgen'), pvar=pvar1)
pvar = fread(paste0(POPS_SNP_data_EGG2, '.pvar'))
psam = fread(paste0(POPS_SNP_data_EGG2, '.psam'))

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

table(is.element(SNPList$dumID1,pvar$ID))
table(is.element(SNPList$dumID2,pvar$ID))

SNPList = SNPList[dumID1 %in% pvar$ID | dumID2 %in% pvar$ID, ]
setorder(SNPList,chr,pos_b38)
table(SNPList$dumID1 == pvar$ID)
table(SNPList$dumID2 == pvar$ID)
table(SNPList$pos_b38 == pvar$POS)
table(SNPList$chr == pvar$CHR)
SNPList[,ID := pvar$ID]

#' # Checks ####
#' ***
#' ## Allele frequency ####
stopifnot(pvar$ID == SNPList$ID)
SNPList[,maf := eaf]
SNPList[eaf>0.5,maf := 1-eaf]

plot(pvar$MAF,SNPList$maf)
plot(pvar$AF,SNPList$eaf)

table(pvar$ALT == SNPList$effect_allele,pvar$REF == SNPList$other_allele)
table(pvar$REF == SNPList$effect_allele,pvar$ALT == SNPList$other_allele)
filt = pvar$ALT == SNPList$effect_allele
plot(pvar$AF[filt],SNPList$eaf[filt])
plot(pvar$AF[!filt],SNPList$eaf[!filt])

#' Okay, there are 8 SNPs with different effect alleles, and one SNP with allele frequency difference. I will switch the alleles (change EAF and beta only), and flag the outlier
SNPList[,EAF2 := eaf]
SNPList[!filt, EAF2 := 1-EAF2]
SNPList[,beta2 := beta]
SNPList[!filt, beta2 := (-1)*beta2]

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
#' How good is the imputation of the selected SNPs?
hist(pvar$R2)
table(pvar$R2<0.8)
pvar[R2<0.8,c(1:5,8:11)]

#' Okay, there are 31 SNPs with low imputation quality. For the moment I will keep them, I can filter for R2 later.  
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

#' Okay, there are no SNPs violating HWE (after correcting for multiple testing). 
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

#' Okay, there are 4271 SNPs with high quality (similar AF as in EGG data, high R2, no sig HWE violation). 
#' 
filt = pvar$highQC == T
pvar = pvar[filt,]
geno_mat = geno_mat[,filt]
SNPList = SNPList[filt,]

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
dev.off()

table(pvar$comment_LD)

#' Okay, there are 17 SNPs with no LD proxy, and 4254 SNP which have at least one LD proxy with LD r2>0.05. 
#' 
#' To reduce the computational burden, I will use some position based pruning, and then check again. 
#' 
myTab = copy(SNPList)
setorder(myTab,p)

getSmallestDist = function(x) {
  if(length(x)>1){
    y = c(x[2:length(x)], max(x)+1000000)
    z = min(y-x)
  }else{
    z=1000000
  }
  return(z)
}


myCHR = unique(myTab$chr)
result.22 = foreach(s2 = myCHR) %do% {
  # s2 = myCHR[1]
  subdata2 = copy(myTab)
  subdata2 = subdata2[chr == s2, ]
  
  setkey(subdata2, pos_b38)
  
  if(dim(subdata2)[1]<=1){
    subdata2[, keep := T]
    subdata2[, NR_SNPs := 0]
  }else{
    subdata2[, keep := NA]
    subdata2[, NR_SNPs := as.numeric(NA)]
    
    smallestDist = getSmallestDist(subdata2[, pos_b38])
    while(smallestDist < 500000) {
      minP = min(subdata2[is.na(keep), p])
      myPOS = subdata2[minP == p & is.na(keep), pos_b38]
      if(length(myPOS)>1){
        myPOS = myPOS[1]
      }
      subdata2[pos_b38 == myPOS, keep := T]
      
      #filter for SNPs that can stay within the set (outside the +- 500 kb range or keep==T)
      myFilt = (subdata2[, pos_b38] < (myPOS - 500000)) | 
        (subdata2[, pos_b38] > (myPOS + 500000)) | 
        subdata2[, keep] 
      myFilt[is.na(myFilt)] = FALSE
      subdata2 = subdata2[myFilt == TRUE, ]
      
      subdata2[pos_b38 == myPOS, NR_SNPs := sum(myFilt==F)]
      smallestDist = getSmallestDist(subdata2[, pos_b38])
    }
    
    #stopifnot(sum(is.na(subdata2[,keep])) <= 1)
    subdata2[is.na(keep), NR_SNPs := 0]
    subdata2[is.na(keep), keep := TRUE]
  }
  
  subdata2
}
SNPList_filtered = rbindlist(result.22)

#' now repeat the LD correlation table
#' 
pvar2 = copy(pvar)
pvar2 = pvar2[ID %in% SNPList_filtered$ID,]
pvar2[,comment_LD := "LD OK"]
filt = is.element(pvar$ID,pvar2$ID)
geno_mat2 = geno_mat[,filt]
myCHRs = unique(pvar2$CHR)

LDTab2 = foreach(i = 1:length(myCHRs))%do%{
  #i=20
  myCHR = myCHRs[i]
  filt = pvar2$CHR == myCHR
  dumMatrix = geno_mat2[,filt]
  if(sum(filt)>1){
    CorTab = cor(dumMatrix)^2
    heatmap(CorTab,Rowv = NA,Colv = NA, main =paste0("Chromosome ",i))
    CorTab2 = as.data.table(CorTab)
    CorTab2[,SNP1 := rownames(CorTab)]
    CorTab_long = melt(CorTab2,id.vars="SNP1",measure.vars=rownames(CorTab),variable.name="SNP2")
    CorTab_long = CorTab_long[SNP1 != SNP2,]
    CorTab2_long = CorTab_long[value>0.05,]
    pvar2[ID %in% CorTab2_long$SNP1 | ID %in% CorTab2_long$SNP2, comment_LD := "LD >0.05 with at least 1 SNP"]
    CorTab_long
  }
  
}
LDTab2 = rbindlist(LDTab2)
dev.off()

table(pvar2$comment_LD)


#' # Save data ####
#' ***
save(pvar, psam, geno_mat, file = paste0(POPS_phenotypes,"/01_Prep_02_SNPData_",tag,".RData"))
save(pvar2, psam, geno_mat2, file = paste0(POPS_phenotypes,"/01_Prep_02_SNPData_filtered_",tag,".RData"))

save(HWETab, file= paste0("../results/01_Prep_02_HWE_",tag,".RData"))
save(LDTab, file = paste0("../results/01_Prep_02_LD_",tag,".RData"))
save(LDTab2, file = paste0("../results/01_Prep_02_LD_filtered_",tag,".RData"))
save(SNPList, file = paste0("../results/01_Prep_02_SNPList_",tag,".RData"))
save(SNPList_filtered, file = paste0("../results/01_Prep_02_SNPList_filtered_",tag,".RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
