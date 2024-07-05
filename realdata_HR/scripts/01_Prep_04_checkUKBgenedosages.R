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
#' In my PLINK2 call, I already filter for imputation quality (mach r2 >0.8) and MAF (>0.01). Now I only have to check the AF again (I do not expect me to filter anything), and if possible to compare to the GWAS Catalog RAF. But RAF is not available for all SNPs, hence I will not filter outliers here.  
#'  
#' Then, I check per chromosome the pairwise LD and perform some crude position-based priority pruning, using the log-transformed p-values from the GWAS Catalog entries for ranking.
#' 
#' This final genedosage information will then be saved for the later GAMLSS runs. 
#' 
#' In addition, I want to load my phenotype data and check the PCs. 
#'  
#' 
#' # Init ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_HPC.R")
.libPaths()

#' # Load genetic data ####
#' ***
UKB_genotypes_filtered = gsub("~","../../../",UKB_genotypes_filtered)
pvar1 = NewPvar(paste0(UKB_genotypes_filtered, 'UKB_HR_merged.pvar'))
pgen = NewPgen(paste0(UKB_genotypes_filtered,'UKB_HR_merged.pgen'), pvar=pvar1)
pvar = fread(paste0(UKB_genotypes_filtered,'UKB_HR_merged.pvar'))
psam = fread(paste0(UKB_genotypes_filtered,'UKB_HR_merged.psam'))

geno_mat <- ReadList(pgen, variant_subset = c(1:dim(pvar)[1]) , meanimpute=F)
dim(geno_mat)
colnames(geno_mat) = pvar$ID
rownames(geno_mat) = psam$IID

head(psam)
setnames(psam, "#FID","FID")

head(pvar)
setnames(pvar, "#CHROM","CHR")
pvar[,EAF := colSums(geno_mat)/(2*30111)]
pvar[,MAF := EAF]
pvar[EAF>0.5,MAF := 1-EAF]
pvar[,table(MAF<0.01)]
save(pvar, psam, geno_mat, file = paste0(UKB_genotypes_filtered,"UKB_HR_merged.RData"))

#' # Some checks ####
#' ***
#' I want to check the following things: 
#' 
#' - EAF similar to EAF reported in the GWAS Catalog?
#' - pairwise LD (which SNPs are actually independent?)
#' 
#' ## AF check #### 
load("../results/01_Prep_03_SNPList_HR.RData")

matched = match(pvar$ID,myTab3$SNPS)
myTab4 = copy(myTab3)
myTab4 = myTab4[matched,]
plot(pvar$EAF,myTab4$RAF)
abline(0,1)
abline(1,-1)
plot(pvar$MAF,myTab4$MAF)
abline(0,1)
abline(1,-1)

#' Okay, there are five outliers, but there are also a lot SNPs with NA as RAF. I still want to flag those weird ones.
#' 
pvar[,comment_AF := "AF OK"]
filt = is.na(myTab4$RAF)
pvar[filt,comment_AF := "AF check not possible (no RAF in GWAS Catalog)"]
filt = abs(myTab4$MAF - pvar$MAF)>0.05
pvar[filt,comment_AF := "AF different"]

#' ## LD check #### 
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

#' Okay, there are 107 SNPs with no LD proxy, and 230 SNP which have at least one LD proxy with LD r2>0.05. 
#' 
#' To reduce the computational burden, I will use some position based pruning, and then check again. 
#' 
pvar[,LOGP_GWASCatalog := myTab4[,PVALUE_MLOG]]

getSmallestDist = function(x) {
  if(length(x)>1){
    y = c(x[2:length(x)], max(x)+1000000)
    z = min(y-x)
  }else{
    z=1000000
  }
  return(z)
}

result.22 = foreach(s2 = myCHRs) %do% {
  # s2 = myCHR[1]
  subdata2 = copy(pvar)
  subdata2 = subdata2[CHR == s2, ]
  setkey(subdata2, POS)
  
  if(dim(subdata2)[1]<=1){
    subdata2[, keep := T]
    subdata2[, NR_SNPs := 0]
  }else{
    subdata2[, keep := NA]
    subdata2[, NR_SNPs := as.numeric(NA)]
    
    smallestDist = getSmallestDist(subdata2[, POS])
    while(smallestDist < 500000) {
      #minP = min(subdata2[is.na(keep), p])
      maxLogP = max(subdata2[is.na(keep), LOGP_GWASCatalog])
      myPOS = subdata2[maxLogP == LOGP_GWASCatalog & is.na(keep), POS]
      if(length(myPOS)>1){
        myPOS = myPOS[1]
      }
      subdata2[POS == myPOS, keep := T]
      
      #filter for SNPs that can stay within the set (outside the +- 500 kb range or keep==T)
      myFilt = (subdata2[, POS] < (myPOS - 500000)) | 
        (subdata2[, POS] > (myPOS + 500000)) | 
        subdata2[, keep] 
      myFilt[is.na(myFilt)] = FALSE
      subdata2 = subdata2[myFilt == TRUE, ]
      
      subdata2[POS == myPOS, NR_SNPs := sum(myFilt==F)]
      smallestDist = getSmallestDist(subdata2[, POS])
    }
    
    #stopifnot(sum(is.na(subdata2[,keep])) <= 1)
    subdata2[is.na(keep), NR_SNPs := 0]
    subdata2[is.na(keep), keep := TRUE]
  }
  
  subdata2
}
SNPList_filtered = rbindlist(result.22)
setorder(SNPList_filtered,CHR,POS)

#' now repeat the LD correlation table
#' 
pvar2 = copy(pvar)
pvar2 = pvar2[ID %in% SNPList_filtered$ID | CHR>20,]
pvar2[,comment_LD := "LD OK"]
filt = is.element(pvar$ID,pvar2$ID)
geno_mat2 = geno_mat[,filt]
test = pvar[,.N,CHR]
test = test[N>1,]
myCHRs = test$CHR

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
save(pvar, psam, geno_mat, file = paste0(UKB_phenotypes_filtered,"/01_Prep_04_SNPData.RData"))
save(pvar2, psam, geno_mat2, file = paste0(UKB_phenotypes_filtered,"/01_Prep_04_SNPData_filtered.RData"))

save(LDTab, file = paste0("../results/01_Prep_04_LD.RData"))
save(LDTab2, file = paste0("../results/01_Prep_04_LD_filtered.RData"))

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
