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
#' Here, I want to run the MVMR using the statistics from all samples with at least 2 exposure measurements. 
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
load("../results/01_Prep_02_LD_240317.RData")
load("../results/01_Prep_02_SNPList_240317.RData")
SNPList[,goodSNP := F]
SNPList[grepl("OK",comment_AF) & grepl("OK",comment_R2) & grepl("OK",comment_HWE),goodSNP := T]
SNPList[,table(goodSNP)]

#' Filter LD table for good SNPs only
LDTab = LDTab[SNP1 %in% SNPList[goodSNP==T,SNP],]
LDTab = LDTab[SNP2 %in% SNPList[goodSNP==T,SNP],]
LDTab[,SNP2 := as.character(SNP2)]

#' ## Exposure
load("../results/02_SNPs_01_MAIN_Assocs_exposure_gamlssIA_all2_gaQuad_240318.RData")
myAssocs_X = copy(myAssocs_X_gamlssIA)
myAssocs_X = myAssocs_X[SNP %in% SNPList[goodSNP==T,SNP]]

#' Check the GX associations
myAssocs_X[,table(pval_mean==0,phenotype)]
myAssocs_X[,table(pval_mean<1e-50,phenotype)]

#' There is some serious issue with the phenotype hc_cm. When revisiting this in the 02_ script, I found convergence problems for this phenotype. None of the other phenotypes is affected. So I remove this one for further analyses. 
#' 
myAssocs_X = myAssocs_X[phenotype != "hc_cm"]

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

#' ## Outcome
load("../results/03_SNPs_01_MAIN_Assocs_outcome_nTIA_240319.RData")
myAssocs_Y = myAssocs_Y[SNP %in% SNPList[goodSNP==T,SNP]]
myAssocs_Y = myAssocs_Y[population == "all",]

#' ## save as temporary files
save(myAssocs_X_long,myAssocs_Y, file = paste0("../temp/04_MVMRInput_MAIN_",tag,".RData"))

#' # Do MVMR ####
#' ***
#' I want 
#' 
#' - all SNPs (including all non-significant ones)
#' - 200 random SNPs (to allow cond. F-statistics estimation)
#' - nominal significant SNPs
#' - (suggestive significant SNPs)
#' - top 20 overlap
#' - top 20 distinct (as good as possible)
#' 
#' I want to do this classically (all three exposures mean, slope and variability), and seperated for mean and variability or slope and variability.
#' 
#' ## Analysis 1: mean, slope and variability
#' 
myExposures = unique(myAssocs_X_long$phenotype)
myOutcomes = unique(myAssocs_Y$phenotype)

registerDoParallel(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))
#registerDoParallel(4)

dumTab2 = foreach(j = 1:length(myExposures))%dorng%{
  #dumTab2 = foreach(j = 1:length(myExposures))%dopar%{
  #j=1
  source("../../SourceFile_HPC.R")
  source("../../helperfunctions/MVMR_jp_POPS.R")
  source("../../helperfunctions/MVMR_jp_POPS_top20.R")
    
  myExposure = myExposures[j]
  
  # filter data
  myAssocs_X_long2 = copy(myAssocs_X_long)
  myAssocs_X_long2 = myAssocs_X_long2[phenotype == myExposure,]
  myAssocs_X_long2[,dumID := "main"]
    
  dumTab3 = foreach(k = 1:length(myOutcomes))%do%{
      #k=1
      myOutcome = myOutcomes[k]
      myAssocs_Y2 = copy(myAssocs_Y)
      myAssocs_Y2 = myAssocs_Y2[phenotype == myOutcome,]
      
      message("Working on exposure ",myExposure," and outcome ",myOutcome," ...")
      
      # do MVMRs
      MVMR0 = MVMR_jp_POPS(data_exposure = myAssocs_X_long2,
                           data_outcome = myAssocs_Y2,
                           exposure_name = myExposure, 
                           outcome_name = myOutcome,
                           flag = "main",
                           GX_pval_treshold = 1,
                           getPlot = F,
                           corTab = LDTab,
                           corTab_threshold = 0.1,sampleSize_GX = 3*2996,random = F,getCondF = F,getUni = T)
      
      MVMR1 = MVMR_jp_POPS(data_exposure = myAssocs_X_long2,
                           data_outcome = myAssocs_Y2,
                           exposure_name = myExposure, 
                           outcome_name = myOutcome,
                           flag = "main",
                           GX_pval_treshold = 1,
                           getPlot = F,
                           corTab = LDTab,
                           corTab_threshold = 0.1,sampleSize_GX = 3*2996,random = T,getCondF = T,getUni = F)
      
      MVMR2 = MVMR_jp_POPS(data_exposure = myAssocs_X_long2,
                           data_outcome = myAssocs_Y2,
                           exposure_name = myExposure, 
                           outcome_name = myOutcome,
                           flag = "main",
                           GX_pval_treshold = 0.05,
                           getPlot = F,
                           corTab = LDTab,
                           corTab_threshold = 0.1,sampleSize_GX = 3*2996,random = F,getCondF = F,getUni = T)
      
      MVMR3 = MVMR_jp_POPS(data_exposure = myAssocs_X_long2,
                           data_outcome = myAssocs_Y2,
                           exposure_name = myExposure, 
                           outcome_name = myOutcome,
                           flag = "main",
                           GX_pval_treshold = 0.05,
                           getPlot = F,
                           corTab = LDTab,
                           corTab_threshold = 0.1,sampleSize_GX = 3*2996,random = T,getCondF = T,getUni = F)
      
      MVMR4 = MVMR_jp_POPS_top20(data_exposure = myAssocs_X_long2,
                                 data_outcome = myAssocs_Y2,
                                 exposure_name = myExposure, 
                                 outcome_name = myOutcome,
                                 flag = "main",
                                 SNPSets = "overlap",
                                 getPlot = F,
                                 corTab = LDTab,
                                 corTab_threshold = 0.1,sampleSize_GX = 3*2996,random = F,getCondF = T,getUni = T)
      
      MVMR5 = MVMR_jp_POPS_top20(data_exposure = myAssocs_X_long2,
                                 data_outcome = myAssocs_Y2,
                                 exposure_name = myExposure, 
                                 outcome_name = myOutcome,
                                 flag = "main",
                                 SNPSets = "distinct",
                                 getPlot = F,
                                 corTab = LDTab,
                                 corTab_threshold = 0.1,sampleSize_GX = 3*2996,random = F,getCondF = T,getUni = T)
      
      MVMR0[,threshold := "no p-value filter, all SNPs"]
      MVMR1[,threshold := "no p-value filter, random 5% of all SNPs"]
      MVMR2[,threshold := "p<0.05, all SNPs"]
      MVMR3[,threshold := "p<0.05, random 5% of all SNPs"]
      MVMR4[,threshold := "top20_overlap"]
      MVMR5[,threshold := "top20_distinct"]
      MVMR = rbind(MVMR0,MVMR1,MVMR2,MVMR3,MVMR4,MVMR5,fill=T)
      MVMR
      
    }
    MVMR_Tab1 = rbindlist(dumTab3)
    MVMR_Tab1
    
}

MVMR_results1 = rbindlist(dumTab2,fill = T)
save(MVMR_results1,file = paste0("../temp/04_MVMR_01_MAIN_classic_",tag,".RData"))

#' ## Analysis 2: mean and variability
#' 
registerDoParallel(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))

dumTab2 = foreach(j = 1:length(myExposures))%dorng%{
  #dumTab2 = foreach(j = 1:length(myExposures))%dopar%{
  #j=6
  source("../../SourceFile_HPC.R")
  source("../../helperfunctions/MVMR_jp_POPS.R")
  source("../../helperfunctions/MVMR_jp_POPS_top20.R")
  
  myExposure = myExposures[j]
  
  # filter data
  myAssocs_X_long2 = copy(myAssocs_X_long)
  myAssocs_X_long2 = myAssocs_X_long2[phenotype == myExposure,]
  myAssocs_X_long2[,dumID := "main"]
  
  dumTab3 = foreach(k = 1:length(myOutcomes))%do%{
    #k=1
    myOutcome = myOutcomes[k]
    myAssocs_Y2 = copy(myAssocs_Y)
    myAssocs_Y2 = myAssocs_Y2[phenotype == myOutcome,]
    
    message("Working on exposure ",myExposure," and outcome ",myOutcome," ...")
    
    # do MVMRs
    MVMR0 = MVMR_jp_POPS(data_exposure = myAssocs_X_long2[type != "slope",],
                         data_outcome = myAssocs_Y2,
                         exposure_name = myExposure, 
                         outcome_name = myOutcome,
                         flag = "main",
                         GX_pval_treshold = 1,
                         getPlot = F,
                         corTab = LDTab,
                         corTab_threshold = 0.1,sampleSize_GX = 3*2996,random = F,getCondF = F,getUni = T)
    
    MVMR1 = MVMR_jp_POPS(data_exposure = myAssocs_X_long2[type != "slope",],
                         data_outcome = myAssocs_Y2,
                         exposure_name = myExposure, 
                         outcome_name = myOutcome,
                         flag = "main",
                         GX_pval_treshold = 1,
                         getPlot = F,
                         corTab = LDTab,
                         corTab_threshold = 0.1,sampleSize_GX = 3*2996,random = T,getCondF = T,getUni = F)
    
    MVMR2 = MVMR_jp_POPS(data_exposure = myAssocs_X_long2[type != "slope",],
                         data_outcome = myAssocs_Y2,
                         exposure_name = myExposure, 
                         outcome_name = myOutcome,
                         flag = "main",
                         GX_pval_treshold = 0.05,
                         getPlot = F,
                         corTab = LDTab,
                         corTab_threshold = 0.1,sampleSize_GX = 3*2996,random = F,getCondF = F,getUni = T)
    
    MVMR3 = MVMR_jp_POPS(data_exposure = myAssocs_X_long2[type != "slope",],
                         data_outcome = myAssocs_Y2,
                         exposure_name = myExposure, 
                         outcome_name = myOutcome,
                         flag = "main",
                         GX_pval_treshold = 0.05,
                         getPlot = F,
                         corTab = LDTab,
                         corTab_threshold = 0.1,sampleSize_GX = 3*2996,random = T,getCondF = T,getUni = F)
    
    MVMR4 = MVMR_jp_POPS_top20(data_exposure = myAssocs_X_long2[type != "slope",],
                               data_outcome = myAssocs_Y2,
                               exposure_name = myExposure, 
                               outcome_name = myOutcome,
                               flag = "main",
                               SNPSets = "overlap",
                               getPlot = F,
                               corTab = LDTab,
                               corTab_threshold = 0.1,sampleSize_GX = 3*2996,random = F,getCondF = T,getUni = T)
    
    MVMR5 = MVMR_jp_POPS_top20(data_exposure = myAssocs_X_long2[type != "slope",],
                               data_outcome = myAssocs_Y2,
                               exposure_name = myExposure, 
                               outcome_name = myOutcome,
                               flag = "main",
                               SNPSets = "distinct",
                               getPlot = F,
                               corTab = LDTab,
                               corTab_threshold = 0.1,sampleSize_GX = 3*2996,random = F,getCondF = T,getUni = T)
    
    MVMR0[,threshold := "no p-value filter, all SNPs"]
    MVMR1[,threshold := "no p-value filter, random 5% of all SNPs"]
    MVMR2[,threshold := "p<0.05, all SNPs"]
    MVMR3[,threshold := "p<0.05, random 5% of all SNPs"]
    MVMR4[,threshold := "top20_overlap"]
    MVMR5[,threshold := "top20_distinct"]
    MVMR = rbind(MVMR0,MVMR1,MVMR2,MVMR3,MVMR4,MVMR5,fill=T)
    MVMR
    
  }
  MVMR_Tab1 = rbindlist(dumTab3)
  MVMR_Tab1
  
}

MVMR_results2 = rbindlist(dumTab2,fill = T)
save(MVMR_results2,file = paste0("../temp/04_MVMR_01_MAIN_mean_",tag,".RData"))

#' ## Analysis 2: mean and variability
#' 
registerDoParallel(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))

dumTab2 = foreach(j = 1:length(myExposures))%dorng%{
  #dumTab2 = foreach(j = 1:length(myExposures))%dopar%{
  #j=6
  source("../../SourceFile_HPC.R")
  source("../../helperfunctions/MVMR_jp_POPS.R")
  source("../../helperfunctions/MVMR_jp_POPS_top20.R")
  
  myExposure = myExposures[j]
  
  # filter data
  myAssocs_X_long2 = copy(myAssocs_X_long)
  myAssocs_X_long2 = myAssocs_X_long2[phenotype == myExposure,]
  myAssocs_X_long2[,dumID := "main"]
  
  dumTab3 = foreach(k = 1:length(myOutcomes))%do%{
    #k=1
    myOutcome = myOutcomes[k]
    myAssocs_Y2 = copy(myAssocs_Y)
    myAssocs_Y2 = myAssocs_Y2[phenotype == myOutcome,]
    
    message("Working on exposure ",myExposure," and outcome ",myOutcome," ...")
    
    # do MVMRs
    MVMR0 = MVMR_jp_POPS(data_exposure = myAssocs_X_long2[type != "mean",],
                         data_outcome = myAssocs_Y2,
                         exposure_name = myExposure, 
                         outcome_name = myOutcome,
                         flag = "main",
                         GX_pval_treshold = 1,
                         getPlot = F,
                         corTab = LDTab,
                         corTab_threshold = 0.1,sampleSize_GX = 3*2996,random = F,getCondF = F,getUni = T)
    
    MVMR1 = MVMR_jp_POPS(data_exposure = myAssocs_X_long2[type != "mean",],
                         data_outcome = myAssocs_Y2,
                         exposure_name = myExposure, 
                         outcome_name = myOutcome,
                         flag = "main",
                         GX_pval_treshold = 1,
                         getPlot = F,
                         corTab = LDTab,
                         corTab_threshold = 0.1,sampleSize_GX = 3*2996,random = T,getCondF = T,getUni = F)
    
    MVMR2 = MVMR_jp_POPS(data_exposure = myAssocs_X_long2[type != "mean",],
                         data_outcome = myAssocs_Y2,
                         exposure_name = myExposure, 
                         outcome_name = myOutcome,
                         flag = "main",
                         GX_pval_treshold = 0.05,
                         getPlot = F,
                         corTab = LDTab,
                         corTab_threshold = 0.1,sampleSize_GX = 3*2996,random = F,getCondF = F,getUni = T)
    
    MVMR3 = MVMR_jp_POPS(data_exposure = myAssocs_X_long2[type != "mean",],
                         data_outcome = myAssocs_Y2,
                         exposure_name = myExposure, 
                         outcome_name = myOutcome,
                         flag = "main",
                         GX_pval_treshold = 0.05,
                         getPlot = F,
                         corTab = LDTab,
                         corTab_threshold = 0.1,sampleSize_GX = 3*2996,random = T,getCondF = T,getUni = F)
    
    MVMR4 = MVMR_jp_POPS_top20(data_exposure = myAssocs_X_long2[type != "mean",],
                               data_outcome = myAssocs_Y2,
                               exposure_name = myExposure, 
                               outcome_name = myOutcome,
                               flag = "main",
                               SNPSets = "overlap",
                               getPlot = F,
                               corTab = LDTab,
                               corTab_threshold = 0.1,sampleSize_GX = 3*2996,random = F,getCondF = T,getUni = T)
    
    MVMR5 = MVMR_jp_POPS_top20(data_exposure = myAssocs_X_long2[type != "mean",],
                               data_outcome = myAssocs_Y2,
                               exposure_name = myExposure, 
                               outcome_name = myOutcome,
                               flag = "main",
                               SNPSets = "distinct",
                               getPlot = F,
                               corTab = LDTab,
                               corTab_threshold = 0.1,sampleSize_GX = 3*2996,random = F,getCondF = T,getUni = T)
    
    MVMR0[,threshold := "no p-value filter, all SNPs"]
    MVMR1[,threshold := "no p-value filter, random 5% of all SNPs"]
    MVMR2[,threshold := "p<0.05, all SNPs"]
    MVMR3[,threshold := "p<0.05, random 5% of all SNPs"]
    MVMR4[,threshold := "top20_overlap"]
    MVMR5[,threshold := "top20_distinct"]
    MVMR = rbind(MVMR0,MVMR1,MVMR2,MVMR3,MVMR4,MVMR5,fill=T)
    MVMR
    
  }
  MVMR_Tab1 = rbindlist(dumTab3)
  MVMR_Tab1
  
}

MVMR_results3 = rbindlist(dumTab2,fill = T)
save(MVMR_results3,file = paste0("../temp/04_MVMR_01_MAIN_slope_",tag,".RData"))

#' # Save data ####
#' ***
MVMR_results = rbind(MVMR_results1,MVMR_results2,MVMR_results3)
save(MVMR_results,file = paste0("../temp/04_MVMR_01_MAIN_",tag,".RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
