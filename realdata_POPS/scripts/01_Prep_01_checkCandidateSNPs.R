#' ---
#' title: "Check candidate SNPs from EGG"
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
#' In the real data application of my longitudinal MVMR analyses, I want to use SNPs affecting (estimated) birth weight.
#' 
#' Here, I load the summary statistics for birth weight from the EGG consortium (downloaded from the [EGG Consortium wegpage](http://egg-consortium.org/birth-weight-2016.html) on 06/02/2024).
#' 
#' "We are releasing the summary data from our genome-wide meta-analyses of birth weight, combining data from the EGG consortium and the first tranche of the UK Biobank data (released May 2015). Data were imputed up to the reference panels from the 1000 Genomes Project (Phase 1 v3) or combined 1000G and UK10K Project. Birth weight was z-score transformed in males and females separately. The association between each variant and birth weight was tested separately in males and females using a linear regression model, with adjustment for gestational week (where available) and study-specific ancestry covariates."
#' 
#' "Summary files provide information on chromosome, genomic position (NCBI build 37), rsID, effect allele, other allele, effect allele frequency, beta, standard error, P value and sample size at over 16 million variants passing quality control."
#' 
#' When using data from the downloadable meta-analyses results please acknowledge the source of the data as follows:
#' 
#' **Data on birth weight trait has been contributed by the EGG Consortium using the UK Biobank Resource and has been downloaded from www.egg-consortium.org.**
#' 
#' In addition to the above acknowledgement, please cite the paper below:
#' 
#' **Genome-wide associations for birth weight and correlations with adult disease. Horikoshi M, Beaumont RN, Day FR, Warrington NM, Kooijman MN, Fernandez-Tajes J, et al. Nature 2016 doi:10.1038/nature19806**
#'  
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile.R")
.libPaths()
suppressPackageStartupMessages(library(SNPlocs.Hsapiens.dbSNP150.GRCh38))

#' # Load EGG data ####
#' ***
EGG_data = fread(EGG_BW, header=TRUE, sep="\t")
head(EGG_data)
names(EGG_data)[2] = "pos_b37"
min(EGG_data$p)
table(EGG_data$p<5e-8)
table(EGG_data$p<1e-6)
table(EGG_data$chr)

#' # Filter EGG data ####
#' ***
#' - pvalue < 1e-6
#' - EAF < 0.99 and EAF > 0.01
#' - autosomal SNPs (chr<23)
#' - rsID available
#' 
EGG_data = EGG_data[p<1e-6,]
EGG_data = EGG_data[eaf>0.01 & eaf<0.99,]
EGG_data = EGG_data[chr<=22,]
EGG_data = EGG_data[grepl("rs",rsid)]

#' # Get hg38 position information ####
#' ***
#' 
mySNPs = EGG_data$rsid
snps = SNPlocs.Hsapiens.dbSNP150.GRCh38
myLift = snpsById(snps, mySNPs, ifnotfound="drop")

myPos = pos(myLift)
myIDs = myLift$RefSNP_id

matched = match(mySNPs,myIDs)
EGG_data[,pos_b38 := myPos[matched]]
EGG_data = EGG_data[!is.na(pos_b38)]

#' # Save data ####
#' ***
#' First, I get the ID in the style required for POPS. Then I create the data file Jasmine can use to extract the SNPs. Finally, I save the EGG data for later check (similar AF in EGG and POPS?)
#' 
EGG_data[,dumID1 := paste0("chr",chr,":",pos_b38,":",effect_allele,":",other_allele)]
EGG_data[,dumID2 := paste0("chr",chr,":",pos_b38,":",other_allele,":",effect_allele)]

outFile1 = paste0(POPS_SNP_data,"01_Prep_01_SNPList_EGG_IDonly.txt")
mySNPlist = c(EGG_data$dumID1,EGG_data$dumID2)
#write.table(mySNPlist, file = outFile1,col.names = F,row.names = F, quote = F)

outFile2 = paste0("../results/01_Prep_01_SNPList.RData")
save(EGG_data, file = outFile2)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
