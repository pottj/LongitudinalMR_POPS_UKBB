#' ---
#' title: "Get MVMR v2"
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
#' Here, I want to run the MVMR using the SNP statistics from the linear mixed model, gamlss, and gamlssIA and test for causal effects of the intercept and slope. 
#' 
#' **I only use the top 20 associations per exposure type!**
#' 
#' Two settings: 
#' 
#' - best 20 SNPs allowing SNP overlaps
#' - best 20 SNPs allowing no SNP overlap
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_HPC.R")
source("../helperfunctions/MVMR_jp_POPS_top20.R")
.libPaths()

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("202.-","24-",tag)
tag = gsub("-","",tag)

#' # Load data ####
#' ***
#' Load Outcome data
load("../results/08_3_Associations_SNPs_outcome_lin_231206.RData")
load("../results/08_7_Associations_SNPs_outcome_glm_240117.RData")

#' Load exposure data
load("../results/08_1_Associations_SNPs_exposure_linMix_231205.RData")
load("../results/08_2_Associations_SNPs_exposure_gamlss_231205.RData")
load("../results/08_5_Associations_SNPs_exposure_gamlss_IA_231211.RData")

#' # Prepare data ####
#' ***
#' Combine outcome data into one data table
myAssocs_Y = rbind(myAsscos_X_linear,myAsscos_X_glm)
myOutcomes = unique(myAssocs_Y$phenotype)

#' Change column names & model in the gamlssIA data table & combine data tables
myAsscos_X_gamlssIA[,model:="gamlssIA"]
names(myAsscos_X_linearMixed)[8:11] = names(myAsscos_X_gamlssIA)[8:11]

myAssocs_X = rbind(myAsscos_X_linearMixed,myAsscos_X_gamlss,myAsscos_X_gamlssIA, fill=T)

#' Transform myAssocs_X into long format
data_long1 = melt(myAssocs_X,
                  id.vars=c("model","SNP","phenotype"),
                  measure.vars=c("beta_mean", "beta_slope", "beta_var"),
                  variable.name="type",
                  value.name="beta")
data_long2 = melt(myAssocs_X,
                  id.vars=c("model","SNP","phenotype"),
                  measure.vars=c("SE_mean", "SE_slope", "SE_var"),
                  variable.name="type",
                  value.name="SE")
data_long3 = melt(myAssocs_X,
                  id.vars=c("model","SNP","phenotype"),
                  measure.vars=c("tval_mean", "tval_slope", "tval_var"),
                  variable.name="type",
                  value.name="tval")
data_long4 = melt(myAssocs_X,
                  id.vars=c("model","SNP","phenotype"),
                  measure.vars=c("pval_mean", "pval_slope", "pval_var"),
                  variable.name="type",
                  value.name="pval")

myAssocs_X_long = cbind(data_long1,data_long2[,5],data_long3[,5],data_long4[,5])
myAssocs_X_long[,type := gsub("beta_","",type)]
myAssocs_X_long = myAssocs_X_long[!is.na(beta),]
setorder(myAssocs_X_long,phenotype,SNP,model,type)

myExposures = unique(myAssocs_X_long$phenotype)

#' ## Test one pair ####
test_overlap = MVMR_jp_POPS_top20(data_exposure=myAssocs_X_long,
                                  data_outcome=myAssocs_Y,
                                  exposure_name=myExposures[11],
                                  outcome_name=myOutcomes[4],
                                  GX_model="gamlssIA",
                                  SNPSets = "overlap",
                                  getPlot=T)

test_distinct = MVMR_jp_POPS_top20(data_exposure=myAssocs_X_long,
                                  data_outcome=myAssocs_Y,
                                  exposure_name=myExposures[11],
                                  outcome_name=myOutcomes[4],
                                  GX_model="gamlssIA",
                                  SNPSets = "distinct",
                                  getPlot=T)

test_overlap
test_distinct

#' Looks good!
#'
#' ## Loop for all MVMRs ###
#' 
#' - Level 1: choose exposure
#' - Level 2: choose outcome
#' - Level 3: choose model
#' 
myModels = unique(myAssocs_X_long$model)

dumTab1 = foreach(i = 1:length(myExposures))%do%{
  # i=1
  myExposure = myExposures[i]
  
  dumTab2 = foreach(j = 1:length(myOutcomes))%do%{
    # j=1
    myOutcome = myOutcomes[j]
    
    dumTab3 = foreach(k = 1:length(myModels))%do%{
      # k=1
      myModel = myModels[k]
      
      test1 = MVMR_jp_POPS_top20(data_exposure=copy(myAssocs_X_long),
                                 data_outcome=copy(myAssocs_Y),
                                 exposure_name=myExposure,
                                 outcome_name=myOutcome,
                                 GX_model=myModel,
                                 SNPSets = "overlap",
                                 getPlot=F)      
      test2 = MVMR_jp_POPS_top20(data_exposure=copy(myAssocs_X_long),
                                 data_outcome=copy(myAssocs_Y),
                                 exposure_name=myExposure,
                                 outcome_name=myOutcome,
                                 GX_model=myModel,
                                 SNPSets = "distinct",
                                 getPlot=F)      
      test1[,threshold := "overlap"]
      test2[,threshold := "distinct"]
      
      test = rbind(test1,test2)
      test
    }
    dumTab3 = rbindlist(dumTab3,fill = T)
    dumTab3
  }
  dumTab2 = rbindlist(dumTab2,fill = T)
  dumTab2
}
MVMR_top20 = rbindlist(dumTab1,fill = T)

save(MVMR_top20, file = "../results/08_8_MVMR.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
