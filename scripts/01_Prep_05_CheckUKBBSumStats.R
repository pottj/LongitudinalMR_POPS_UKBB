#' ---
#' title: "Prepare UKBB data"
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
#' I want to run the MVMR in a real 2-sample MR setting using POPS for the exposure effects, and UKBB for the outcome effects (self-reported birth weight).
#' 
#' I downloaded the data from [Neale lab](http://www.nealelab.is/uk-biobank), and here I want to prepare the data for later use. I need to load the summary statistics and match the SNP position and alleles. POPS uses GRCh38, while the UKBB is in GRCh37. However, the PGS data has both the harmonized GRCh38 and the original GRCh37 positions. So I load them all and filter accordingly. I will check if the effect alleles are the same, and how similar the effect allele frequencies are. Finally, I will save the data for later use. 
#' 
#' # Initialization ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_HPC.R")

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("202.-","24-",tag)
tag = gsub("-","",tag)

#' # Load data ####
#' ***
#' ## Load POPS ####
#' GRCh38 coordinates
pvar = fread("../../rds/rds-obsgynae-POPs-DC0Q0rbedlk/POPs_Imputation/Pott/Pott_Fetal_PGS.pvar")
head(pvar)

#' Split the information in the last column
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

head(pvar)
pvar[,INFO:=NULL]

#' ## Load PGS data ####
#' GRCh37 and GRCh38 coordinates
load("../temp/06_PGSInputData.RData")
PGS_data = matchedData
table(is.element(pvar$ID,PGS_data$SNP))
PGS_data = PGS_data[SNP %in% pvar$ID]
stopifnot(pvar$ID == PGS_data$SNP)

#' ## Load UKBB SumStats ####
#' GRCh37 coordinates
#' 
#' I only load the summary statistics for the *raw* phenotype data here, the data for the inverse-rank normal transformed phenotype, *irnt*, will be loaded later. The reason for this is that they have the same structure and variant ID, so I only have to do this once, and can later simply use the *raw* data table to match the *irnt* one. 
#'   
GWAS_raw=read.table("../../data/downloadedData/20022_raw.gwas.imputed_v3.both_sexes.tsv.bgz",header=T,sep="\t") 
setDT(GWAS_raw)
GWAS_raw = GWAS_raw[!is.na(pval)]
GWAS_raw = GWAS_raw[low_confidence_variant == "false"]
dim(GWAS_raw)
head(GWAS_raw)

#' # Match UKB and PGS ####
#' ***
#' The previous dumIDs were for *hm_pos*, which gives the base position in GRCh38. Now I need other dumIDs with the base positions in GRCh37, as given in *chr_position*.
PGS_data[,table(hm_pos == chr_position)]
PGS_data[,dumID3 := paste0(chr_name,":",chr_position,":",effect_allele,":",other_allele)]
PGS_data[,dumID4 := paste0(chr_name,":",chr_position,":",other_allele,":",effect_allele)]

table(is.element(PGS_data$dumID3,GWAS_raw$variant))
table(is.element(PGS_data$dumID4,GWAS_raw$variant))

#' Okay, only *dumID4* is relevant, which makes sense as the readme from Neale lab said that the second allele is the effect allele in their GWAS (so effect alleles should already be a match!)
#' 
GWAS_raw = GWAS_raw[variant %in% PGS_data$dumID4,]
PGS_data = PGS_data[dumID4 %in% GWAS_raw$variant,]
pvar = pvar[ID %in% PGS_data$SNP,]

#' # Checks ####
#' ***
#' ## Check 1: same order ####
table(GWAS_raw$variant == PGS_data$dumID4)
matched = match(PGS_data$dumID4,GWAS_raw$variant)
GWAS_raw = GWAS_raw[matched,]
table(GWAS_raw$variant == PGS_data$dumID4)

#' ## Check 2: same effect allele ####
#' I want additional columns for effect allele and other allele
#' 
dummy = GWAS_raw$variant
dummy2 = unlist(strsplit(dummy,":"))

dummy_chr = dummy2[seq(from=1,to=23264,by=4)]
stopifnot(PGS_data$chr_name == dummy_chr)
dummy_pos = dummy2[seq(from=2,to=23264,by=4)]
stopifnot(PGS_data$chr_position == dummy_pos)
dummy_OA = dummy2[seq(from=3,to=23264,by=4)]
dummy_EA = dummy2[seq(from=4,to=23264,by=4)]

GWAS_raw[,effect_allele := dummy_EA]
GWAS_raw[,other_allele := dummy_OA]

#' Check if the alleles are the same (should be, as matched via dumID4, but best to double check here)
#' 
table(GWAS_raw$effect_allele == pvar$ALT)
table(GWAS_raw$other_allele == pvar$REF)

table(GWAS_raw$effect_allele == PGS_data$`ALT(1)`)
table(GWAS_raw$other_allele == PGS_data$`REF(0)`)

#' ## Check 3: same effect allele frequency ####
#' Now get effect allele frequency
GWAS_raw[effect_allele == minor_allele,EAF := minor_AF]
GWAS_raw[effect_allele != minor_allele,EAF := 1-minor_AF]

plot(GWAS_raw$EAF,pvar$AF)
plot(GWAS_raw$minor_AF,pvar$MAF)
plot(GWAS_raw$EAF,PGS_data$ALT_Frq)
plot(GWAS_raw$minor_AF,PGS_data$MAF)

#' There are some SNPs with weird discrepancy. I remove them to be sure of the matching alleles, as I do not know where the error occured. 
#' 
filt1 = GWAS_raw$EAF<0.5 & pvar$AF>0.5
filt2 = GWAS_raw$EAF>0.5 & pvar$AF<0.5
table(filt1,filt2)

GWAS_raw = GWAS_raw[!filt1 & !filt2,]
pvar = pvar[!filt1 & !filt2,]
PGS_data = PGS_data[!filt1 & !filt2,]

plot(GWAS_raw$EAF,pvar$AF)
plot(GWAS_raw$EAF,PGS_data$ALT_Frq)

#' OK, now it looks good
#' 
#' # Repeat for *irnt* data ####
#' 
GWAS_irnt=read.table("../../data/downloadedData/20022_irnt.gwas.imputed_v3.both_sexes.tsv.bgz",header=T,sep="\t") 
setDT(GWAS_irnt)
GWAS_irnt = GWAS_irnt[variant %in% GWAS_raw$variant,]
matched = match(GWAS_raw$variant,GWAS_irnt$variant)
GWAS_irnt = GWAS_irnt[matched,]
stopifnot(GWAS_irnt$variant == GWAS_raw$variant)

#' Add allele information
GWAS_irnt[,effect_allele := GWAS_raw$effect_allele]
GWAS_irnt[,other_allele := GWAS_raw$other_allele]
GWAS_irnt[effect_allele == minor_allele,EAF := minor_AF]
GWAS_irnt[effect_allele != minor_allele,EAF := 1-minor_AF]

#' Make some plots checking allele frequencies
plot(GWAS_irnt$EAF,pvar$AF)
plot(GWAS_irnt$minor_AF,pvar$MAF)
plot(GWAS_raw$EAF,GWAS_irnt$EAF)
plot(GWAS_raw$minor_AF,GWAS_irnt$minor_AF)

#' Make some plots comparing the beta effect estimates 
plot(GWAS_irnt$beta,GWAS_raw$beta)
plot(GWAS_irnt$beta,PGS_data$effect_weight)
plot(GWAS_raw$beta,PGS_data$effect_weight)

#' # Save data ####
#' ***
#'
GWAS_raw[,SNPID_POPS := pvar$ID]
save(GWAS_raw,file = "../data/UKBB_GWAS_BW_raw.RData")

GWAS_irnt[,SNPID_POPS := pvar$ID]
save(GWAS_irnt,file = "../data/UKBB_GWAS_BW_irnt.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
