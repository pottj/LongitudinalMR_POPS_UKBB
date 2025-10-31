#' ---
#' title: "Create FUMA input"
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
#' I want to annotate the nominal significant SNPs using FUMA. I want per SNP the minimal p-value across all settings (main and sensitivity tests)
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile.R")

#' # Load and prep data ####
#' ***
load("../results/_tables/SupplementalTables_EGG.RData")

myTab = copy(tab2)
myTab[,minPval := pmin(POPS_EFW_pval_mean,POPS_EFW_pval_slope,POPS_EFW_pval_var,na.rm = T)]
myTab[,minPvalType := "mean"]
myTab[minPval == POPS_EFW_pval_slope, minPvalType := "slope"]
myTab[minPval == POPS_EFW_pval_var, minPvalType := "variability"]
myTab[,table(minPvalType)]

setorder(myTab,minPval)
myTab = myTab[!duplicated(rsID)]
myTab[,table(minPvalType,flag)]
myTab[,table(minPvalType,POPS_EFW_phenotype)]

setorder(myTab,chr,pos_b38)

#' FUMA needs only the following columns: 
#' 
#' - chromosome
#' - position (GRCh38 or GRCh37)
#' - rsID
#' - p-value
#' - effect allele
#' - non-effect allele
#' - OR or beta
#' - SE
#' - Sample size
#' 
myTab[minPvalType=="mean",beta := POPS_EFW_beta_mean]
myTab[minPvalType=="slope",beta := POPS_EFW_beta_slope]
myTab[minPvalType=="variability",beta := POPS_EFW_beta_var]

myTab[minPvalType=="mean",SE := POPS_EFW_SE_mean]
myTab[minPvalType=="slope",SE := POPS_EFW_SE_slope]
myTab[minPvalType=="variability",SE := POPS_EFW_SE_var]

names(myTab)
myNames = names(myTab)[c(3,4,2,40,6,7,42,43,10)]
myNames
FUMA = copy(myTab)
colsOut<-setdiff(colnames(FUMA),myNames)
FUMA[,get("colsOut"):=NULL]
setcolorder(FUMA,myNames)
names(FUMA) = c("chr","pos","rsID","pval","EA","OA","beta","SE","N")

#' # Save ####
#' ***
fwrite(FUMA,file = "../results/_tables/FUMA_input.txt",sep = "\t")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
