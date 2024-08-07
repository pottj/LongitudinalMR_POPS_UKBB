#' ---
#' title: "MVMR - main analysis"
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
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_HPC.R")
source("../../helperfunctions/MVMR_jp_POPS.R")
source("../../helperfunctions/MVMR_jp_POPS_top20.R")

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("202.-","24-",tag)
tag = gsub("-","",tag)

#' # Get data ####
#' ***
#' ## SNP info
#' Get high quality SNPs only
load("../results/01_Prep_04_LD_filtered.RData")
LDTab = copy(LDTab2)
LDTab[,SNP2 := as.character(SNP2)]
load("../results/01_Prep_03_SNPList_HR.RData")

#' ## Exposure
load("../results/02_SNPs_01_MAIN.RData")

#' Check the GX associations
myAssocs_X[,table(pval_mean==0)]
myAssocs_X[,table(pval_mean<1e-50)]
myTab3[SNPS %in% myAssocs_X[pval_mean<1e-50,SNP]]
myAssocs_X[,table(pval_slope<5e-8)]
myTab3[SNPS %in% myAssocs_X[pval_slope<5e-8,SNP]]
myAssocs_X[,table(pval_var<1e-30)]
myTab3[SNPS %in% myAssocs_X[pval_var<1e-30,SNP]]

myAssocs_X[,table(pval_mean<5e-8,pval_slope<5e-8)]
cor.test(myAssocs_X$beta_mean,myAssocs_X$beta_slope)
cor.test(myAssocs_X$beta_mean,myAssocs_X$beta_var)
cor.test(myAssocs_X$beta_slope,myAssocs_X$beta_var)

#' Transform into wide format
data_long1 = melt(myAssocs_X,
                  id.vars=names(myAssocs_X)[1:3],
                  measure.vars=c("beta_mean", "beta_slope", "beta_var"),
                  variable.name="type",
                  value.name="beta")
data_long2 = melt(myAssocs_X,
                  id.vars=names(myAssocs_X)[1:3],
                  measure.vars=c("SE_mean", "SE_slope", "SE_var"),
                  variable.name="type",
                  value.name="SE")
data_long3 = melt(myAssocs_X,
                  id.vars=names(myAssocs_X)[1:3],
                  measure.vars=c("tval_mean", "tval_slope", "tval_var"),
                  variable.name="type",
                  value.name="tval")
data_long4 = melt(myAssocs_X,
                  id.vars=names(myAssocs_X)[1:3],
                  measure.vars=c("pval_mean", "pval_slope", "pval_var"),
                  variable.name="type",
                  value.name="pval")

myAssocs_X_long = cbind(data_long1,data_long2[,5],
                        data_long3[,5],data_long4[,5])
myAssocs_X_long[,type := gsub("beta_","",type)]
myAssocs_X_long = myAssocs_X_long[!is.na(beta),]

#' ## Outcome
load("../results/01_Prep_06_CADsummaryStats.RData")

#' ## save as temporary files
goodSNPs = myAssocs_Y[,.N,by = rsID]
goodSNPs = goodSNPs[N==4,]
myAssocs_X_long = myAssocs_X_long[SNP %in% goodSNPs$rsID,]
myAssocs_Y = myAssocs_Y[rsID %in% goodSNPs$rsID,]
save(myAssocs_X_long,myAssocs_Y, file = paste0("../temp/03_MVMRInput_MAIN.RData"))

#' # Do MVMR ####
#' ***
#' I want 
#' 
#' - all SNPs (including all non-significant ones)
#' - nominal significant SNPs
#' - top 20 overlap
#' - top 20 distinct (as good as possible)
#' 
#' I want to do this classically (all three exposures mean, slope and variability), and seperated for mean and variability or slope and variability.
#' 
#' ## Analysis 1: mean, slope and variability
#' 
myExposures = unique(myAssocs_X_long$model)
mySampleSize = c(22349, 11309, 11040)
myOutcomes = unique(myAssocs_Y$phenotype)
myFlag = "main"

setnames(myAssocs_Y,"SNP","markername")
setnames(myAssocs_Y,"rsID","SNP")

registerDoParallel(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))
#registerDoParallel(4)

#dumTab2 = foreach(j = 1:length(myExposures))%dorng%{
dumTab2 = foreach(j = 1:length(myExposures))%dopar%{
  #j=2
  source("../../SourceFile_HPC.R")
  source("../../helperfunctions/MVMR_jp_POPS.R")
  source("../../helperfunctions/MVMR_jp_POPS_top20.R")
    
  myExposure = myExposures[j]
  
  # filter data
  myAssocs_X_long2 = copy(myAssocs_X_long)
  myAssocs_X_long2 = myAssocs_X_long2[model == myExposure,]
  myAssocs_X_long2[,phenotype := paste0("HR_",myExposure)]
  myAssocs_X_long2[,dumID := myFlag]
    
  dumTab3 = foreach(k = 1:length(myOutcomes))%do%{
      #k=2
      myOutcome = myOutcomes[k]
      myAssocs_Y2 = copy(myAssocs_Y)
      myAssocs_Y2 = myAssocs_Y2[phenotype == myOutcome,]
      
      message("Working on exposure ",myExposure," and outcome ",myOutcome," ...")
      
      # do MVMRs
      MVMR0 = MVMR_jp_POPS(data_exposure = myAssocs_X_long2,
                           data_outcome = myAssocs_Y2,
                           exposure_name = paste0("HR_",myExposure), 
                           outcome_name = myOutcome,
                           flag = myFlag,
                           GX_pval_treshold = 1,
                           getPlot = F,
                           corTab = LDTab,
                           corTab_threshold = 0.1,sampleSize_GX = mySampleSize[j],
                           random = F,getCondF = T,getUni = T)

      MVMR2 = MVMR_jp_POPS(data_exposure = myAssocs_X_long2,
                           data_outcome = myAssocs_Y2,
                           exposure_name = paste0("HR_",myExposure), 
                           outcome_name = myOutcome,
                           flag = myFlag,
                           GX_pval_treshold = 1e-6,
                           getPlot = F,
                           corTab = LDTab,
                           corTab_threshold = 0.1,sampleSize_GX = mySampleSize[j],
                           random = F,getCondF = T,getUni = T)

      MVMR4 = MVMR_jp_POPS_top20(data_exposure = myAssocs_X_long2,
                                 data_outcome = myAssocs_Y2,
                                 exposure_name = paste0("HR_",myExposure), 
                                 outcome_name = myOutcome,
                                 flag = myFlag,
                                 SNPSets = "overlap",
                                 getPlot = F,
                                 corTab = LDTab,
                                 corTab_threshold = 0.1,sampleSize_GX = mySampleSize[j],
                                 random = F,getCondF = T,getUni = T)
      
      MVMR5 = MVMR_jp_POPS_top20(data_exposure = myAssocs_X_long2,
                                 data_outcome = myAssocs_Y2,
                                 exposure_name = paste0("HR_",myExposure), 
                                 outcome_name = myOutcome,
                                 flag = myFlag,
                                 SNPSets = "distinct",
                                 getPlot = F,
                                 corTab = LDTab,
                                 corTab_threshold = 0.1,sampleSize_GX = mySampleSize[j],
                                 random = F,getCondF = T,getUni = T)
      
      MVMR0[,threshold := "all_SNPs"]
      MVMR2[,threshold := "nominal_SNPs"]
      MVMR4[,threshold := "top20_overlap"]
      MVMR5[,threshold := "top20_distinct"]
      MVMR = rbind(MVMR0,MVMR2,MVMR4,MVMR5,fill=T)
      MVMR
      
    }
    MVMR_Tab1 = rbindlist(dumTab3)
    MVMR_Tab1
    
}

MVMR_results = rbindlist(dumTab2,fill = T)
save(MVMR_results,file = paste0("../results/03_MVMR_01_MAIN.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
