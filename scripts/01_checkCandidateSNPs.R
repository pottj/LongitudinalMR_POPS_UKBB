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
data = fread("../data/gwas-association-downloaded_2023-09-22-EFO_0004344.tsv")

#' Reduce to relevant columns
#' 
names(data)
data = data[,c(2:4,9:15,22,21,27,31,28)]
names(data) = c("PMID","Author1","Date","SampleSize_init","SampleSize_rep","cytoband","chr","pos","genes_reported","genes_mapped","rsID","EA","EAF","BETA","pval")
data[,EA := gsub(".*[-]","",EA)]

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

#' Check for allele frequency
data[,min(EAF)]
data = data[EAF>=0.005,]
data = data[EAF<=0.995,]

#' Check remaining cytobands
#' 
test = data[,.N, by=cytoband]
table(test$N)

#' There are 205 unique entries, of which 107 only occur once. There are some loci with multiple SNPs - here I keep all those SNPs, and select later the best-associated in the POPs analysis. 
#' 
#' # Get SNP ID from UK Biobank ### 
#' ***
#' Jasmine asked for ID in style of CHR#:POS:REF:ALT, e.g. chr3:159837169:T:G. 
#' 
#' REF = Effect allele (also sometimes called risk allele, reference allele, effect allele, coded allele, etc.)
#' 
#' ALT = Non-effect allele (also some times called alternate allele, the other allele etc.)
#' 
#' Here, I use as REF the risk allele as given in the GWAS catalog data. The ALT allele is the not-given allele as listed in the dbSNP database. I remove all tri-allelic SNPs to make my life easier. 
#'  
library(SNPlocs.Hsapiens.dbSNP150.GRCh38)
snps = SNPlocs.Hsapiens.dbSNP150.GRCh38
my_rsids = data$rsID
my_snps = snpsById(snps, my_rsids)
dummy1 = my_snps@elementMetadata@listData$alleles_as_ambig
dummy2 = my_snps@elementMetadata@listData$RefSNP_id
dummy = data.table::data.table(ID = dummy2, coding = dummy1)

#' - R = A or G	
#' - K = G or T	
#' - S = G or C
#' - Y = C or T	
#' - M = A or C	
#' - W = A or T

dummy[coding == "R", a1 := "A"]
dummy[coding == "R", a2 := "G"]

dummy[coding == "K", a1 := "G"]
dummy[coding == "K", a2 := "T"]

dummy[coding == "S", a1 := "G"]
dummy[coding == "S", a2 := "C"]

dummy[coding == "Y", a1 := "C"]
dummy[coding == "Y", a2 := "T"]

dummy[coding == "M", a1 := "A"]
dummy[coding == "M", a2 := "C"]

dummy[coding == "W", a1 := "A"]
dummy[coding == "W", a2 := "T"]

dummy = dummy[!is.na(a2)]
data = data[rsID %in% dummy$ID,]
stopifnot(data$rsID == dummy$ID)
table(data$EA == dummy$a1, data$EA == dummy$a2)

data[,OA := dummy$a2]
data[,OA2 := dummy$a1]
data[OA == EA, OA := OA2]
data[,OA2 := NULL]
table(data$EA == data$OA)

data[, chr_pos_REF_ALT := paste0("chr",chr,":",pos,":",EA,":",OA)]
data[, chr_pos_ALT_REF := paste0("chr",chr,":",pos,":",OA,":",EA)]
data = data[,c(1:12,16,17,18,13:15)]

#' # Save ####
#' ***
outFile = paste0("../results/01_SNPList_BirthWeight_",tag,".txt")
fwrite(data, file = outFile)

outFile2 = paste0("../results/01_SNPList_BirthWeight_SNPList_",tag,".txt")
list = c(data$chr_pos_REF_ALT,data$chr_pos_ALT_REF)
table(duplicated(list))
write.table(list, file = outFile2,col.names = F,row.names = F, quote = F)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
