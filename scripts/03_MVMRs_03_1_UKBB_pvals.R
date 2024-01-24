#' ---
#' title: "Get MVMR with UKBB outcome data"
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
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_HPC.R")
source("../helperfunctions/MVMR_jp_POPS.R")
.libPaths()

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("202.-","24-",tag)
tag = gsub("-","",tag)

#' # Load data ####
#' ***
#' Load Outcome data
load("../data/UKBB_GWAS_BW_raw.RData")
load("../data/UKBB_GWAS_BW_irnt.RData")

#' Load exposure data
load("../results/08_1_Associations_SNPs_exposure_linMix_231205.RData")
load("../results/08_2_Associations_SNPs_exposure_gamlss_231205.RData")
load("../results/08_5_Associations_SNPs_exposure_gamlss_IA_231211.RData")

#' # Prepare data ####
#' ***
#' Combine outcome data into one data table
GWAS_irnt[,phenotype := "UKBB_BW_irnt"]
GWAS_raw[,phenotype := "UKBB_BW_raw"]
myAssocs_Y = rbind(GWAS_raw,GWAS_irnt)
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

#' ## Filter SNPs ####
#' In the UKBB, I could not match every SNP. Hence I have to filter my exposure data for those SNPs only. 
#' 
myAssocs_X = myAssocs_X[SNP %in% myAssocs_Y$SNPID_POPS,]
myAssocs_X_long = myAssocs_X_long[SNP %in% myAssocs_Y$SNPID_POPS,]
names(myAssocs_X)
names(myAssocs_Y)

setnames(myAssocs_Y,"SNPID_POPS","SNP")
setnames(myAssocs_Y,"beta","beta_mean")
setnames(myAssocs_Y,"se","SE_mean")
setnames(myAssocs_Y,"tstat","tval_mean")
setnames(myAssocs_Y,"pval","pval_mean")

#' ## Test one pair ####
test = MVMR_jp_POPS(data_exposure=myAssocs_X_long,
                    data_outcome=myAssocs_Y,
                    exposure_name=myExposures[11],
                    outcome_name=myOutcomes[1],
                    GX_model="gamlssIA",
                    GX_pval_treshold = 1e-3,
                    getPlot=T)
test

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
  # i=4
  myExposure = myExposures[i]
  
  dumTab2 = foreach(j = 1:length(myOutcomes))%do%{
    # j=1
    myOutcome = myOutcomes[j]
    
    dumTab3 = foreach(k = 1:length(myModels))%do%{
      # k=3
      myModel = myModels[k]
      
      test1 = MVMR_jp_POPS(data_exposure=myAssocs_X_long,
                           data_outcome=myAssocs_Y,
                           exposure_name=myExposure,
                           outcome_name=myOutcome,
                           GX_model=myModel,
                           GX_pval_treshold = 0.05,
                           getPlot=F)
      
      test2 = MVMR_jp_POPS(data_exposure=myAssocs_X_long,
                           data_outcome=myAssocs_Y,
                           exposure_name=myExposure,
                           outcome_name=myOutcome,
                           GX_model=myModel,
                           GX_pval_treshold = 1e-3,
                           getPlot=F)
      
      test1[,threshold := 0.05]
      test2[,threshold := 1e-3]
      
      test = rbind(test1,test2,fill=T)
      test
    }
    dumTab3 = rbindlist(dumTab3,fill = T)
    dumTab3
  }
  dumTab2 = rbindlist(dumTab2,fill = T)
  dumTab2
}
MVMR_UKBB_pvals = rbindlist(dumTab1,fill = T)

save(MVMR_UKBB_pvals, file = "../results/11_3_MVMR_UKBB_pvals.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
