#' ---
#' title: "Check candidate SNPs"
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
#' Here, I load the list of known GWAS hits for birth weight as listed in the GWAS Catalog (EFO 0004344, downloaded on 22/09/2023). 
#' 
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

#' # Load GWAS Catalog data ####
#' ***
#' Complete download of all associations listed for birth weight (EFO ID: EFO_0004344) on Friday, 22/09/2023. 
#' 
data = fread("../temp/gwas-association-downloaded_2023-09-22-EFO_0004344.tsv")

#' Reduce to relevant columns
#' 
names(data)
data = data[,c(2:4,9:15,22,27,31,28)]
names(data) = c("PMID","Author1","Date","SampleSize_init","SampleSize_rep","cytoband","chr","pos","genes_reported","genes_mapped","rsID","EAF","BETA","pval")

#' Reformat some columns
#'
table(data$SampleSize_init)
data[,SampleSize_init := gsub("up to ","",SampleSize_init)]
data[,SampleSize_init := gsub("Up to ","",SampleSize_init)]
data[,SampleSize_init := gsub(" .*","",SampleSize_init)]
data[,SampleSize_init := gsub(",","",SampleSize_init)]
data[,SampleSize_init := as.numeric(SampleSize_init)]

table(data$SampleSize_rep)
data[,SampleSize_rep := gsub("up to ","",SampleSize_rep)]
data[,SampleSize_rep := gsub("Up to ","",SampleSize_rep)]
data[,SampleSize_rep := gsub(" .*","",SampleSize_rep)]
data[,SampleSize_rep := gsub(",","",SampleSize_rep)]
data[,SampleSize_rep := as.numeric(SampleSize_rep)]

table(data$chr)
data[chr=="X",chr:="23"]
data[,chr := as.numeric(chr)]

data[,EAF:= as.numeric(EAF)]
hist(data$EAF)

#' Filter for complete SNP information (chr, pos, EAF, rsID)
#' 
data = data[!is.na(chr) & !is.na(pos) & !is.na(EAF) & !is.na(rsID), ]

#' Check for significant association (p<5e-8)
#' 
data[,max(pval)]
data = data[pval<5e-8,]

#' Check for duplicated SNP IDs (keep entry with larger sample size)
data[, table(duplicated(rsID))]
setorder(data,-SampleSize_init)
data = data[!duplicated(rsID)]

#' Check remaining cytobands
#' 
test = data[,.N, by=cytoband]
table(test$N)

#' There are 206 unique entries, of which 107 only occur once. There are some loci with multiple SNPs - here I keep all those SNPs, and select later the best-associated in the POPs analysis. 
#' 
#' # Save ####
#' ***
outFile = paste0("../results/01_SNPList_BirthWeight_",tag,".txt")
fwrite(data, file = outFile)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
