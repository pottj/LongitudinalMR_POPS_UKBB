#' ---
#' title: "SNP selection real data 1"
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
#' 1) Load data from some consortium (EGG) --> n0
#' 2) Filter SNPs                          --> n1
#'    - p<5x10^-8
#'    - MAF>=1%
#'    - autosomal
#'    - rsID available
#'    - match to hg38 possible
#' 3) Extract genotypes from study (POPS)  --> n2
#' 4) Extract summary statistics from other study (UKB) and restrict to overlap --> n3
#' 5) Priority pruning based on p-values from consortium (EGG) and position (1 MB window) --> n4 = n_final
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_HPC.R")
.libPaths()
suppressPackageStartupMessages(library(SNPlocs.Hsapiens.dbSNP150.GRCh38))
source("../../helperfunctions/getSmallestDist.R")

#' # Load EGG data ####
#' ***
EGG_data = fread(paste0(pathData,"/BW3_EUR_summary_stats.txt.gz"), header=TRUE, sep="\t")
n0 = dim(EGG_data)[1]
head(EGG_data)
names(EGG_data)[2] = "pos_b37"
min(EGG_data$p)
table(EGG_data$p<5e-8)
table(EGG_data$chr)

#' # Filter EGG data ####
#' ***
#' - pvalue < 5e-8
#' - EAF < 0.99 and EAF > 0.01
#' - autosomal SNPs (chr<23)
#' - rsID available
#' - match to hg38 possible
#' 
EGG_data = EGG_data[p<5e-8,]
EGG_data = EGG_data[eaf>0.01 & eaf<0.99,]
EGG_data = EGG_data[chr<=22,]
EGG_data = EGG_data[grepl("rs",rsid)]
EGG_data[,dumID := paste(chr,pos_b37,sep=":")]
dupPos = EGG_data[duplicated(dumID)]
EGG_data = EGG_data[!is.element(dumID,dupPos$dumID)]
EGG_data[,dumID := NULL]

mySNPs = EGG_data$rsid
snps = SNPlocs.Hsapiens.dbSNP150.GRCh38
myLift = snpsById(snps, mySNPs, ifnotfound="drop")

myPos = pos(myLift)
myIDs = myLift$RefSNP_id

matched = match(mySNPs,myIDs)
EGG_data[,pos_b38 := myPos[matched]]
EGG_data = EGG_data[!is.na(pos_b38)]

n1 = dim(EGG_data)[1]

#' Add ID as in POPS data
EGG_data[,dumID1 := paste0("chr",chr,":",pos_b38,":",effect_allele,":",other_allele)]
EGG_data[,dumID2 := paste0("chr",chr,":",pos_b38,":",other_allele,":",effect_allele)]

#' Add ID as in Neale lab data
EGG_data[,dumID3 := paste0(chr,":",pos_b37,":",effect_allele,":",other_allele)]
EGG_data[,dumID4 := paste0(chr,":",pos_b37,":",other_allele,":",effect_allele)]

#' # Extract SNPs in exposure data ####
#' ***
POPS_SNP_data_EGG2 = gsub("~","../../..",POPS_SNP_data_EGG)
pvar1 = NewPvar(paste0(POPS_SNP_data_EGG2, '.pvar'))
pgen = NewPgen(paste0(POPS_SNP_data_EGG2, '.pgen'), pvar=pvar1)
pvar = fread(paste0(POPS_SNP_data_EGG2, '.pvar'))
psam = fread(paste0(POPS_SNP_data_EGG2, '.psam'))

pvar = pvar[ID %in% c(EGG_data$dumID1,EGG_data$dumID2),]
EGG_data = EGG_data[dumID1 %in% pvar$ID | dumID2 %in% pvar$ID,]
n2 = dim(EGG_data)[1]

#' # Extract SNPs in outcome data ####
#' ***
UKB_raw = fread(paste0(pathData,"20022_raw.gwas.imputed_v3.both_sexes.tsv.bgz"))
UKB_raw = UKB_raw[low_confidence_variant ==F,]
table(is.element(EGG_data$dumID3,UKB_raw$variant))
table(is.element(EGG_data$dumID4,UKB_raw$variant))

UKB_raw = UKB_raw[variant %in% EGG_data$dumID3 | variant %in% EGG_data$dumID4,]
EGG_data = EGG_data[dumID3 %in% UKB_raw$variant | dumID4 %in% UKB_raw$variant,]
n3 = dim(EGG_data)[1]

#' # Perform priority pruning ####
#' ***
setorder(EGG_data,chr,pos_b37)
myCHRs = unique(EGG_data$chr)
EGG_data[,pvalue_neg_log10 := -log10(p)]

result.22 = foreach(s2 = myCHRs) %do% {
  # s2 = myCHRs[1]
  subdata2 = copy(EGG_data)
  subdata2 = subdata2[chr == s2, ]
  setkey(subdata2, pos_b37)
  
  if(dim(subdata2)[1]<=1){
    subdata2[, keep := T]
    subdata2[, NR_SNPs := 0]
  }else{
    subdata2[, keep := NA]
    subdata2[, NR_SNPs := as.numeric(NA)]
    
    smallestDist = getSmallestDist(subdata2[, pos_b37])
    while(smallestDist < 1000000) {
      #minP = min(subdata2[is.na(keep), p])
      maxLogP = max(subdata2[is.na(keep), pvalue_neg_log10])
      myPOS = subdata2[maxLogP == pvalue_neg_log10 & is.na(keep), pos_b37]
      if(length(myPOS)>1){
        myPOS = myPOS[1]
      }
      subdata2[pos_b37 == myPOS, keep := T]
      
      #filter for SNPs that can stay within the set (outside the +- 1MB range or keep==T)
      myFilt = (subdata2[, pos_b37] < (myPOS - 1000000)) | 
        (subdata2[, pos_b37] > (myPOS + 1000000)) | 
        subdata2[, keep] 
      myFilt[is.na(myFilt)] = FALSE
      subdata2 = subdata2[myFilt == TRUE, ]
      
      subdata2[pos_b37 == myPOS, NR_SNPs := sum(myFilt==F)]
      smallestDist = getSmallestDist(subdata2[, pos_b37])
    }
    
    #stopifnot(sum(is.na(subdata2[,keep])) <= 1)
    subdata2[is.na(keep), NR_SNPs := 0]
    subdata2[is.na(keep), keep := TRUE]
  }
  
  subdata2
}
EGG_data2 = rbindlist(result.22)
n4 = dim(EGG_data2)[1]

#' # Compare Allele frequencies ####
#' ***
#' ## Checking POPS data
pvar = pvar[ID %in% c(EGG_data2$dumID1,EGG_data2$dumID2),]
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

table(EGG_data2$dumID1==pvar$ID,EGG_data2$dumID2==pvar$ID)
table(EGG_data2$effect_allele==pvar$ALT, EGG_data2$other_allele==pvar$REF)
plot(EGG_data2$eaf, pvar$AF)
abline(0,1)

#' ## Checking UKB data
UKB_raw = UKB_raw[variant %in% EGG_data2$dumID3 | variant %in% EGG_data2$dumID4,]
UKB_raw[,EA := gsub(".*:","",variant)]
UKB_raw[,OA := gsub("[1234567890]","",variant)]
UKB_raw[,OA := gsub("::","",OA)]
UKB_raw[,OA := gsub(":.*","",OA)]
UKB_raw[,EAF := minor_AF]
UKB_raw[EA != minor_allele,EAF := 1-minor_AF]

table(EGG_data2$dumID3==UKB_raw$variant,EGG_data2$dumID4==UKB_raw$variant)
table(EGG_data2$effect_allele==UKB_raw$EA,EGG_data2$other_allele==UKB_raw$OA)
plot(EGG_data2$eaf, UKB_raw$EAF)
abline(0,1)
plot(EGG_data2$beta, UKB_raw$beta)
abline(0,1)
cor.test(EGG_data2$beta, UKB_raw$beta)
plot(EGG_data2$beta/EGG_data2$se, UKB_raw$tstat)
abline(0,1)

#' # Combine and save data 
#' ***
#' I want some finale file including 
#' 
#' - SNP information (rsid, chr, pos b37, pos b38, effect allele, other allele)
#' - matching ID for POPS and POPS eaf
#' - EGG summary statistics (eaf, beta, se, t, p, n)
#' - UKB summary statistics (eaf, beta, se, t, p, n)
#' 
names(EGG_data2)
names(pvar)
names(UKB_raw)

SNPList = data.table(rsID = EGG_data2$rsid,
                     chr = EGG_data2$chr, 
                     pos_b37 = EGG_data2$pos_b37,
                     pos_b38 = EGG_data2$pos_b38, 
                     effect_allele = EGG_data2$effect_allele, 
                     other_allele = EGG_data2$other_allele)

SNPList[, POPS_ID := pvar$ID]
SNPList[, POPS_EAF := pvar$AF]

SNPList[, EGG_EAF := EGG_data2$eaf]
SNPList[, EGG_beta := EGG_data2$beta]
SNPList[, EGG_SE := EGG_data2$se]
SNPList[, EGG_tval := EGG_data2$beta / EGG_data2$se]
SNPList[, EGG_pval := EGG_data2$p]
SNPList[, EGG_sampleSize := EGG_data2$n]

SNPList[, UKB_ID := UKB_raw$variant]
SNPList[, UKB_EAF := UKB_raw$EAF]
SNPList[, UKB_beta := UKB_raw$beta]
SNPList[, UKB_SE := UKB_raw$se]
SNPList[, UKB_tval := UKB_raw$tstat]
SNPList[, UKB_pval := UKB_raw$pval]
SNPList[, UKB_sampleSize := UKB_raw$n_complete_samples]

save(SNPList,file = "../results/01_Prep_05_SNPList.RData")

n0;n1;n2;n3;n4

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
