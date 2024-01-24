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
#' Here, I want to run the MVMR using the SNP statistics from the linear mixed model and gamlss and test for causal effects of the intercept and slope. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_HPC.R")
.libPaths()

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("202.-","23-",tag)
tag = gsub("-","",tag)

#' # Load data ####
#' ***
load("../results/08_1_Associations_SNPs_exposure_linMix_231205.RData")
load("../results/08_2_Associations_SNPs_exposure_gamlss_231205.RData")
load("../results/08_3_Associations_SNPs_outcome_lin_231206.RData")
load("../results/08_5_Associations_SNPs_exposure_gamlss_IA_231211.RData")

#' # Do MVMR ####
#' ***
myExposures = unique(myAsscos_X_linearMixed$phenotype)
myOutcomes = unique(myAsscos_X_linear$phenotype)

myAsscos_X_linearMixed[,table(phenotype,pval_mean<0.05)]
myAsscos_X_linearMixed[,table(phenotype,pval_var<0.05)]
myThresholds_lmm = c(0.05,0.01,0.005)

myAsscos_X_gamlss[,table(phenotype,pval_mean<0.05)]
myAsscos_X_gamlss[,table(phenotype,pval_var<0.05)]
myThresholds_gam = c(0.05,1e-4,5e-6)

myAsscos_X_gamlssIA[,table(phenotype,pval_mean<0.05)]
myAsscos_X_gamlssIA[,table(phenotype,pval_var<0.05)]
myAsscos_X_gamlssIA[,table(phenotype,pval_slope<0.05)]

#' ## Define function ####
MVMR_jp_POPS = function(data_exposure,exposure_name, data_outcome, outcome_name, GX_pval_treshold,doPlotting=F,NR_exposure_types){
  #debug
  # data_exposure = copy(myAsscos_X_gamlss)
  # data_outcome = copy(myAsscos_X_linear)
  # exposure_name = myExposures[2]
  # outcome_name = myOutcomes[1]
  # GX_pval_treshold = myThresholds_gam[2]
  # NR_exposure_types = 2
  
  data_GX_wide = copy(data_exposure)
  data_GX_wide = data_GX_wide[phenotype == exposure_name,]
  if(NR_exposure_types==2){
    data_GX_wide = data_GX_wide[pval_var<GX_pval_treshold | pval_mean<GX_pval_treshold,]
  }else{
    data_GX_wide = data_GX_wide[pval_var<GX_pval_treshold | pval_mean<GX_pval_treshold | pval_slope<GX_pval_treshold,]
    }
  data_GY = copy(data_outcome)
  data_GY = data_GY[phenotype == outcome_name,]
  data_GY = data_GY[SNP %in% data_GX_wide$SNP,]
  
  filt1 = grepl("beta",names(data_GX_wide))
  data_beta = copy(data_GX_wide)
  data_beta = data_beta[,filt1,with=F]
  filt2 = grepl("SE",names(data_GX_wide))
  data_SE = copy(data_GX_wide)
  data_SE = data_SE[,filt2,with=F]
  
  if(dim(data_beta)[1]<3){
    res3 = data.table(exposure = rep(exposure_name,2),
                      exposure_type = c("mean","var"),
                      outcome = rep(outcome_name,2),
                      NR_SNPs = rep(dim(data_GY)[1],2))
  }else{
    if(unique(data_GX_wide$model)=="gamlss") exposure_type_names = c("mean","var")
    if(unique(data_GX_wide$model)=="linearMixed") exposure_type_names = c("mean","slope")
    if(NR_exposure_types==3) exposure_type_names = c("mean","slope","var")
    mvmr_obj = mr_mvinput(bx = as.matrix(data_beta),
                          bxse = as.matrix(data_SE),
                          by = data_GY$beta_mean, 
                          byse = data_GY$SE_mean,
                          exposure = exposure_type_names,
                          outcome = outcome_name)
    
    res2 = mr_mvivw(mvmr_obj)    
    res3 = data.table(exposure = rep(exposure_name,NR_exposure_types),
                      exposure_type = c(res2@Exposure),
                      outcome = rep(res2@Outcome,NR_exposure_types),
                      NR_SNPs = rep(dim(data_GY)[1],NR_exposure_types),
                      beta_IVW = c(res2@Estimate),
                      SE_IVW = c(res2@StdError),
                      pval_IVW = c(res2@Pvalue),
                      HeteroStat = rep(res2@Heter.Stat[1],NR_exposure_types),
                      HeteroStat_pval = rep(res2@Heter.Stat[2],NR_exposure_types))
    
    if(doPlotting==T){
      # make a plot
      plotData = data.table(myID = c(data_GX_wide$SNP,data_GX_wide$SNP),
                            myX = c(data_GX_wide$beta_mean,data_GX_wide$beta_var),
                            myY = rep(data_GY$beta_mean,2),
                            myP_X = c(data_GX_wide$pval_mean,data_GX_wide$pval_var),
                            myP_Y = rep(data_GY$pval_mean,2),
                            exposure = rep(c("mean","var"),each=dim(data_GX_wide)[1]))
      SNPs1 = plotData[myP_X<GX_pval_treshold & exposure == "mean",myID]
      SNPs2 = plotData[myP_X<GX_pval_treshold & exposure != "mean",myID]
      plotData[myID %in% SNPs1,color := 1]
      plotData[myID %in% SNPs2,color := 2]
      plotData[myID %in% SNPs1 & myID %in% SNPs2,color := 3]
      table(plotData$color)
      
      data_hline = data.table(exposure = c("mean","var"),
                              line = c(res3$beta_IVW1))
      data_hline  
      
      plot = ggplot(data = plotData, aes(x = myX, y = myY, color = as.factor(color))) + 
        geom_point() +
        facet_wrap(~exposure,scales = "free")+        
        geom_abline(data = data_hline, aes(slope = line,intercept=0)) +
        scale_colour_manual(values=c("#B2182B","#2166AC","#82B446"),
                            labels=c("mean","var","both")) +
        labs(x=paste0("SNP effect on ",exposure_name), 
             y=paste0("SNP effect on ",outcome_name),
             color="Significant in ")
      
      print(plot)
      
    }
    
  }
  
  return(res3)
  
  
}

#' ## Test one pair ####
test = MVMR_jp_POPS(data_exposure = myAsscos_X_gamlss, 
                    data_outcome = myAsscos_X_linear,
                    exposure_name = myExposures[2],
                    outcome_name = myOutcomes[1],
                    GX_pval_treshold = myThresholds_gam[2],
                    NR_exposure_types = 2)

test

#' ## Define loop for lmm ###
#' 
#' - Level 1: choose exposure
#' - Level 2: choose outcome
#' - Level 3: choose pvalue threshold
#' 

dumTab1 = foreach(i = 1:length(myExposures))%do%{
  # i=1
  myExposure = myExposures[i]
  
  dumTab2 = foreach(j = 1:length(myOutcomes))%do%{
    # j=1
    myOutcome = myOutcomes[j]
    
    dumTab3 = foreach(k = 1:length(myThresholds_lmm))%do%{
      # k=1
      myThreshold = myThresholds_lmm[k]
      
      test = MVMR_jp_POPS(data_exposure = myAsscos_X_linearMixed, 
                          data_outcome = myAsscos_X_linear,
                          exposure_name = myExposure,
                          outcome_name = myOutcome,
                          GX_pval_treshold = myThreshold,
                          NR_exposure_types = 2)
      
      test[,threshold := myThreshold]
      test
    }
    dumTab3 = rbindlist(dumTab3,fill = T)
    dumTab3
  }
  dumTab2 = rbindlist(dumTab2,fill = T)
  dumTab2
}
MVMR_lmm = rbindlist(dumTab1,fill = T)

#' ## Define loop for gam ###
#' 
#' - Level 1: choose exposure
#' - Level 2: choose outcome
#' - Level 3: choose pvalue threshold
#' 

dumTab1 = foreach(i = 1:length(myExposures))%do%{
  # i=2
  myExposure = myExposures[i]
  
  dumTab2 = foreach(j = 1:length(myOutcomes))%do%{
    # j=1
    myOutcome = myOutcomes[j]
    
    dumTab3 = foreach(k = 1:length(myThresholds_gam))%do%{
      # k=2
      myThreshold = myThresholds_gam[k]
      
      test = MVMR_jp_POPS(data_exposure = myAsscos_X_gamlss, 
                          data_outcome = myAsscos_X_linear,
                          exposure_name = myExposure,
                          outcome_name = myOutcome,
                          GX_pval_treshold = myThreshold,
                          NR_exposure_types = 2)
      
      test[,threshold := myThreshold]
      test
    }
    dumTab3 = rbindlist(dumTab3,fill = T)
    dumTab3
  }
  dumTab2 = rbindlist(dumTab2,fill = T)
  dumTab2
}
MVMR_gam = rbindlist(dumTab1,fill = T)

#' ## Define loop for gam ###
#' 
#' - Level 1: choose exposure
#' - Level 2: choose outcome
#' - Level 3: choose pvalue threshold
#' 

dumTab1 = foreach(i = 1:length(myExposures))%do%{
  # i=1
  myExposure = myExposures[i]
  
  dumTab2 = foreach(j = 1:length(myOutcomes))%do%{
    # j=1
    myOutcome = myOutcomes[j]
    
    dumTab3 = foreach(k = 1:length(myThresholds_gam))%do%{
      # k=1
      myThreshold = myThresholds_gam[k]
      
      test = MVMR_jp_POPS(data_exposure = myAsscos_X_gamlssIA, 
                          data_outcome = myAsscos_X_linear,
                          exposure_name = myExposure,
                          outcome_name = myOutcome,
                          GX_pval_treshold = myThreshold,
                          NR_exposure_types = 3)
      
      test[,threshold := myThreshold]
      test
    }
    dumTab3 = rbindlist(dumTab3,fill = T)
    dumTab3
  }
  dumTab2 = rbindlist(dumTab2,fill = T)
  dumTab2
}
MVMR_gamIA = rbindlist(dumTab1,fill = T)

#' ## Combine and save
#'
MVMR_gam[,method:="gamlss"]
MVMR_gamIA[,method:="gamlssIA"]
MVMR_lmm[,method:="linearMixed"]

MVMR = rbind(MVMR_lmm,MVMR_gam,MVMR_gamIA)
save(MVMR, file = "../results/08_6_MVMR.RData")

#' # Forest Plots ####
#' ***
#' I want forest plots for the primary exposures (3), on the primary outcomes (2). 
library(grid)
library(forestploter)

myExposures = myExposures[c(7,12,13)]
myExposures2 = c("EFW (absolute)","EFW (Z-score)","EFW (Centile)")
myOutcomes = myOutcomes[c(4,5)]
myOutcomes2 = c("BW (Centile)","BW (absolute)")
myMethods = unique(MVMR$method)

dumTab4 = foreach(i=1:length(myExposures))%do%{
  #i=1
  myExposure = myExposures[i]
  myExposure2 = myExposures2[i]
    
  data1 = copy(MVMR)
  data1 = data1[exposure == myExposure]
  setorder(data1,exposure_type, outcome)
  data1 = data1[!is.na(beta_IVW)]
  
  dumTab5 = foreach(j=1:length(myOutcomes))%do%{
    #j=1
    myOutcome = myOutcomes[j]
    myOutcome2 = myOutcomes2[j]
    
    data2 = copy(data1)
    data2 = data2[outcome == myOutcome,]
    
    dumTab6 = foreach(k=1:length(myMethods))%do%{
      #k=1
      myMethod = myMethods[k]
      
      data3 = copy(data2)
      data3 = data3[method == myMethod,]
      
      myExposureTypes = unique(data3$exposure_type)
      
      dumTab7 = foreach(l=1:length(myExposureTypes))%do%{
        #l=1
        myExposureType = myExposureTypes[l]
        
        data4 = copy(data3)
        data4 = data4[exposure_type == myExposureType,]
        data4[,subgroup := paste0("GX p-value treshold: ", threshold)]
        data4[,lowerCI95 := beta_IVW-1.96*SE_IVW]
        data4[,upperCI95 := beta_IVW+1.96*SE_IVW]
        data4$` ` <- paste(rep(" ", 50), collapse = " ")
        data4$`Estimate [95% CI]` <- ifelse(is.na(data4$SE_IVW), "",
                                            sprintf("%.2f [%.2f, %.2f]",
                                                    data4$beta_IVW, data4$lowerCI95, data4$upperCI95))
        setnames(data4,"subgroup",paste(myMethod,myExposureType, sep=" - "))
        
        min_data4 = floor(min(data4$lowerCI95))
        max_data4 = ceiling(max(data4$upperCI95))
        
        range = max_data4 - min_data4
        myTicks = seq(min_data4,max_data4,by = range/4)
        
        myXlab = paste0(myExposure2," on ",myOutcome2)
        
        p2<- forest(data4[,c(12,15,16)],
                    est = data4$beta_IVW,
                    lower = data4$lowerCI95, 
                    upper = data4$upperCI95,
                    sizes = 0.5,
                    ticks_at = myTicks,
                    ci_column = 2,
                    ref_line = 0,
                    xlim = c(min_data4, max_data4),
                    xlab = myXlab)
        
        #plot(p2)
        
        filename = paste0("../figures/08_6_ForestPlots_",myExposure,"_",myOutcome,"_",myMethod,"_",myExposureType,"_",tag,".tiff")
        tiff(filename = filename,width = 1700, height = 400, res=200, compression = 'lzw')
        plot(p2)
        dev.off()
        print(myXlab)
      }
    }
  }
}

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
