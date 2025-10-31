#' ---
#' title: "Get MVMR within POPS - main analysis"
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
#' Here, I want to run the MVMR using the statistics from all samples with at least 2 exposure measurements. I want to run only
#' 
#' - exposure = log-transformed EFW
#' - outcomes = un-transformed BW
#'    - in POPS adjusted for GA
#'    - in UKB "raw"
#' - using the 52 selected instruments
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile.R")
source("../../helperfunctions/MVMR_jp_POPS.R")

#' # Get data ####
#' ***
#' ## SNP info
#' Get high quality SNPs only
load("../results/01_Prep_02_LD_filtered.RData")
LDTab = copy(LDTab2)
load("../results/01_Prep_05_SNPList.RData")

#' Filter LD table for good SNPs only
LDTab[,SNP2 := as.character(SNP2)]

#' ## Exposure
load("../results/02_SNPs_01_MAIN.RData")
myAssocs_X = copy(myAssocs_X_gamlssIA)
myAssocs_X = myAssocs_X[rsID %in% SNPList$rsID,]

#' Check the GX associations
myAssocs_X[,cor.test(beta_mean,beta_slope),by=phenotype]
myAssocs_X[,cor.test(beta_mean,beta_var),by=phenotype]
myAssocs_X[,cor.test(beta_slope,beta_var),by=phenotype]

#' Transform into wide format
data_long1 = melt(myAssocs_X,
                  id.vars=names(myAssocs_X)[1:7],
                  measure.vars=c("beta_mean", "beta_slope", "beta_var"),
                  variable.name="type",
                  value.name="beta")
data_long2 = melt(myAssocs_X,
                  id.vars=names(myAssocs_X)[1:7],
                  measure.vars=c("SE_mean", "SE_slope", "SE_var"),
                  variable.name="type",
                  value.name="SE")
data_long3 = melt(myAssocs_X,
                  id.vars=names(myAssocs_X)[1:7],
                  measure.vars=c("tval_mean", "tval_slope", "tval_var"),
                  variable.name="type",
                  value.name="tval")
data_long4 = melt(myAssocs_X,
                  id.vars=names(myAssocs_X)[1:7],
                  measure.vars=c("pval_mean", "pval_slope", "pval_var"),
                  variable.name="type",
                  value.name="pval")

myAssocs_X_long = cbind(data_long1,data_long2[,9],
                        data_long3[,9],data_long4[,9])
myAssocs_X_long[,type := gsub("beta_","",type)]
myAssocs_X_long = myAssocs_X_long[!is.na(beta),]

#' Correct the slope estimates for average age in POPS data: 40.316 weeks
myAssocs_X_long[type=="slope", beta := beta * 40.316]
myAssocs_X_long[type=="slope", SE := SE * 40.316]

#' ## Outcome
load("../results/03_SNPs_01_MAIN_Assocs_outcome_nTIA.RData")
myAssocs_Y = myAssocs_Y[population == "all",]
myAssocs_Y = myAssocs_Y[rsID %in% SNPList$rsID,]
myAssocs_Y = myAssocs_Y[phenotype == "pn_bw",]
myAssocs_Y[,phenotype := "POPS_BW"]

myAssocs_Y2 = copy(myAssocs_Y)
stopifnot(myAssocs_Y2$rsID == SNPList$rsID)
myAssocs_Y2[, phenotype := "UKB_BW"]
myAssocs_Y2[, sampleSize := SNPList$UKB_sampleSize]
myAssocs_Y2[, beta_mean := SNPList$UKB_beta]
myAssocs_Y2[, SE_mean := SNPList$UKB_SE]
myAssocs_Y2[, tval_mean := SNPList$UKB_tval]
myAssocs_Y2[, pval_mean := SNPList$UKB_pval]

myAssocs_Y3 = copy(myAssocs_Y)
stopifnot(myAssocs_Y3$rsID == SNPList$rsID)
myAssocs_Y3[, phenotype := "EGG_BW"]
myAssocs_Y3[, sampleSize := SNPList$EGG_sampleSize]
myAssocs_Y3[, beta_mean := SNPList$EGG_beta]
myAssocs_Y3[, SE_mean := SNPList$EGG_SE]
myAssocs_Y3[, tval_mean := SNPList$EGG_tval]
myAssocs_Y3[, pval_mean := SNPList$EGG_pval]

myAssocs_Y = rbind(myAssocs_Y,myAssocs_Y2,myAssocs_Y3)

#' ## save as temporary files
save(myAssocs_X_long,myAssocs_Y, file = paste0("../temp/04_MVMRInput_MAIN.RData"))

#' # Do MVMR ####
#' ***
#' I want 
#' 
#' - all SNPs (including all non-significant ones)
#' - nominal significant SNPs
#' 
myExposures = unique(myAssocs_X_long$phenotype)
myOutcomes = unique(myAssocs_Y$phenotype)
myFlag = "main"
myAssocs_X_long[,dumID := myFlag]

dumTab2 = foreach(j = 1:length(myExposures))%do%{
  #j=1
  myExposure = myExposures[j]
  myAssocs_X_long2 = copy(myAssocs_X_long)
  myAssocs_X_long2 = myAssocs_X_long2[phenotype == myExposure,]
  
  dumTab3 = foreach(k = 1:length(myOutcomes))%do%{
    #k=1
    myOutcome = myOutcomes[k]
    myAssocs_Y2 = copy(myAssocs_Y)
    myAssocs_Y2 = myAssocs_Y2[phenotype == myOutcome,]
    
    # do MVMRs
    MVMR0 = MVMR_jp_POPS(data_exposure = myAssocs_X_long2,
                         data_outcome = myAssocs_Y2,
                         exposure_name = myExposure, 
                         outcome_name = myOutcome,
                         flag = myFlag,
                         GX_pval_treshold = 1,
                         getPlot = F,
                         corTab = LDTab,
                         corTab_threshold = 0.1,sampleSize_GX = 2996,
                         random = F,getCondF = T,getUni = T)
    
    MVMR2 = MVMR_jp_POPS(data_exposure = myAssocs_X_long2,
                         data_outcome = myAssocs_Y2,
                         exposure_name = myExposure, 
                         outcome_name = myOutcome,
                         flag = myFlag,
                         GX_pval_treshold = 0.05,
                         getPlot = F,
                         corTab = LDTab,
                         corTab_threshold = 0.1,sampleSize_GX = 2996,
                         random = F,getCondF = T,getUni = T)
    
    MVMR0[,threshold := "all_SNPs"]
    MVMR2[,threshold := "nominal_SNPs"]
    MVMR = rbind(MVMR0,MVMR2,fill=T)
    MVMR
    
  }
  MVMR_results2 = rbindlist(dumTab3,fill = T)
  MVMR_results2
}
MVMR_results = rbindlist(dumTab2,fill = T)
save(MVMR_results,file = paste0("../results/04_MVMR_01_MAIN.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
