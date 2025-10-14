#' ---
#' title: "Create FUMA input"
#' subtitle: "Longitudinal MVMR in UKB"
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
#' I want to annotate the genome-wide significant SNPs using FUMA. I want per SNP the minimal p-value across all settings (main and sensitivity tests)
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile.R")

#' # Load and prep data ####
#' ***
load("../results/_tables/SupplementalTables_GLGC.RData")

myTab = copy(tab2)
myTab[,minPval := pmin(UKB_TC_pval_mean,UKB_TC_pval_slope,UKB_TC_pval_var,na.rm = T)]
myTab[,minPvalType := "mean"]
myTab[minPval == UKB_TC_pval_slope, minPvalType := "slope"]
myTab[minPval == UKB_TC_pval_var, minPvalType := "variability"]
myTab[,table(minPvalType)]

setorder(myTab,minPval)
myTab = myTab[!duplicated(rsID)]
myTab[,table(minPvalType,flag)]

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
myTab[minPvalType=="mean",beta := UKB_TC_beta_mean]
myTab[minPvalType=="slope",beta := UKB_TC_beta_slope]
myTab[minPvalType=="variability",beta := UKB_TC_beta_var]

myTab[minPvalType=="mean",SE := UKB_TC_SE_mean]
myTab[minPvalType=="slope",SE := UKB_TC_SE_slope]
myTab[minPvalType=="variability",SE := UKB_TC_SE_var]

names(myTab)
myNames = names(myTab)[c(3,4,2,45,6,7,47,48,10)]
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
