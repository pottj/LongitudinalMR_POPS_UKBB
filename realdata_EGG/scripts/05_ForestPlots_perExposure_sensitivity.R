#' ---
#' title: "Get Sensitivity Forest Plots"
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
#' **Forest Plots with sensitivity data**
#' 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../../SourceFile_HPC.R")

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("2024-","24-",tag)
tag = gsub("-","",tag)

#' # Get data ####
#' ***
matchingTable = data.table(old = c("efwcomb","logefwcomb","efwcombZv2","efwcombv2_cent","pn_bw","BW_Centile_Br1990","BW_SDS_Br1990","pn_emcsall"),
                           new = c("EFW_R","EFW_L","EFW_Z","EFW_C","BW_R","BW_C","BW_Z","eCS"))

load("../results/04_MVMR_01_MAIN.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
matched = match(MVMR_results$exposure,matchingTable$old)
MVMR_results[,exposure := matchingTable[matched,new]]
matched = match(MVMR_results$outcome,matchingTable$old)
MVMR_results[,outcome := matchingTable[matched,new]]
MVMR_main = copy(MVMR_results)

load("../results/04_MVMR_02_SENS_GBR3.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
matched = match(MVMR_results$exposure,matchingTable$old)
MVMR_results[,exposure := matchingTable[matched,new]]
matched = match(MVMR_results$outcome,matchingTable$old)
MVMR_results[,outcome := matchingTable[matched,new]]
MVMR_GBR3 = copy(MVMR_results)
MVMR_GBR3[,ID := gsub("sens_","sens1_",ID)]

load("../results/04_MVMR_03_SENS_noVar.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
matched = match(MVMR_results$exposure,matchingTable$old)
MVMR_results[,exposure := matchingTable[matched,new]]
matched = match(MVMR_results$outcome,matchingTable$old)
MVMR_results[,outcome := matchingTable[matched,new]]
MVMR_noVar = copy(MVMR_results)
MVMR_noVar[,ID := gsub("sens_","sens2_",ID)]

load("../results/04_MVMR_04_SENS_noSlope.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
matched = match(MVMR_results$exposure,matchingTable$old)
MVMR_results[,exposure := matchingTable[matched,new]]
matched = match(MVMR_results$outcome,matchingTable$old)
MVMR_results[,outcome := matchingTable[matched,new]]
MVMR_noSlope = copy(MVMR_results)
MVMR_noSlope[,ID := gsub("sens_","sens3_",ID)]

load("../results/04_MVMR_05_SENS_SigmaTimeIndep.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
matched = match(MVMR_results$exposure,matchingTable$old)
MVMR_results[,exposure := matchingTable[matched,new]]
matched = match(MVMR_results$outcome,matchingTable$old)
MVMR_results[,outcome := matchingTable[matched,new]]
MVMR_SigmaTimeIndep = copy(MVMR_results)
MVMR_SigmaTimeIndep[,ID := gsub("sens_","sens4_",ID)]

load("../results/04_MVMR_06_SENS_noCovars.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
matched = match(MVMR_results$exposure,matchingTable$old)
MVMR_results[,exposure := matchingTable[matched,new]]
matched = match(MVMR_results$outcome,matchingTable$old)
MVMR_results[,outcome := matchingTable[matched,new]]
MVMR_noCovars = copy(MVMR_results)
MVMR_noCovars[,ID := gsub("sens_","sens5_",ID)]

load("../results/04_MVMR_07_SENS_it200.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
matched = match(MVMR_results$exposure,matchingTable$old)
MVMR_results[,exposure := matchingTable[matched,new]]
matched = match(MVMR_results$outcome,matchingTable$old)
MVMR_results[,outcome := matchingTable[matched,new]]
MVMR_it200 = copy(MVMR_results)
MVMR_it200[,ID := gsub("sens_","sens6_",ID)]

MVMR = rbind(MVMR_main,MVMR_GBR3,MVMR_noVar,MVMR_noSlope,MVMR_SigmaTimeIndep,MVMR_noCovars,MVMR_it200)
MVMR[,table(ID)]

#' # Create Forest Plots ####
#' ***
#' I need to do this separately for each phenotype and each exposure type. Let's try a loop... 
#' 
myOutcomes = unique(MVMR$outcome)
myExposures = unique(MVMR$exposure)
myExposure_Types = unique(MVMR$exposure_type)

dumTab0 = foreach(k = 1:length(myOutcomes))%do%{
  #k=1
  myOutcome = myOutcomes[k]
  message("Working on outcome ",myOutcome)
  MVMR0 = copy(MVMR)
  MVMR0 = MVMR0[outcome == myOutcome,]
  
  dumTab1 = foreach(i = 1:length(myExposures))%do%{
    #i=1
    myExposure = myExposures[i]
    message("   and exposure ",myExposure)
    MVMR1 = copy(MVMR0)
    MVMR1 = MVMR1[exposure == myExposure]
    
    dumTab2 = foreach(j = 1:length(myExposure_Types))%do%{
      #j=1
      myExposure_Type = myExposure_Types[j]
      message("         and type ",myExposure_Type)
      
      MVMR2 = copy(MVMR1)
      MVMR2 = MVMR2[exposure_type == myExposure_Type]
      MVMR2 = MVMR2[!grepl("random",threshold),]
      
      # start plotting 
      MVMR2[,rank := 2]
      dummy = data.table(ID = unique(MVMR2$ID),rank=1)
      data4 = rbind(MVMR2,dummy, fill=T)
      data4[,subgroup := paste0("   ", threshold)]
      data4[is.na(threshold),subgroup := ID]
      data4[,subgroup := gsub("_", " - ",subgroup)]
      setnames(data4,"subgroup", "   Instrument selection")
      
      data4[,lowerCI95 := beta_IVW-1.96*SE_IVW]
      data4[,upperCI95 := beta_IVW+1.96*SE_IVW]
      data4$` ` <- paste(rep(" ", 50), collapse = " ")
      data4$`Estimate [95% CI]` <- ifelse(is.na(data4$SE_IVW), "",
                                          sprintf("%.2f [%.2f, %.2f]",
                                                  data4$beta_IVW, data4$lowerCI95, data4$upperCI95))
      
      min_data4 = round(min(data4$lowerCI95,na.rm = T) - 0.1,2)
      max_data4 = round(max(data4$upperCI95,na.rm = T) + 0.1,2)
      range = max_data4 - min_data4
      myTicks = seq(min_data4, max_data4, by = round(range/5,2))
      
      setorder(data4,ID,rank)
      # data4 = data4[c(15:21,8:14,1:7,22:28),]
      #data4[threshold == "nominal",threshold := "p<0.05, all SNPs"]
      
      dummy = data4$threshold
      dummy[is.na(dummy)] = "white"
      dummy[grepl("all_SNPs",dummy)] = "lightgrey"
      dummy[grepl("nominal_SNPs",dummy)] = "lightsteelblue2"
      dummy = gsub("top20_overlap","lightgreen",dummy)
      dummy = gsub("top20_distinct","lightcoral",dummy)
      tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))
      
      myXlab = paste0(myExposure_Type, " of ",myExposure," on ",myOutcome)
      
      p2<- forest(data4[,c(16,19,20)],
                  est = data4$beta_IVW,
                  lower = data4$lowerCI95, 
                  upper = data4$upperCI95,
                  sizes = 0.5,
                  ticks_at = myTicks,ticks_digits = 2,
                  #ticks_at = c(0,1,2,3),ticks_digits = 2,
                  ci_column = 2,
                  ref_line = 0,
                  xlim = c(min_data4, max_data4),
                  #xlim = c(-0.5, 3.5),
                  title = myXlab,
                  theme = tm1)
      
      plot(p2)
      
      myX = gsub("_","",myExposure)
      myY = gsub("_","",myOutcome)
      filename = paste0("../results/_figures/05_ForestPlots_sens/",myX,"_",myY,"_",myExposure_Type,".png")
      png(filename = filename,width = 2000, height = 2000, res=200)
      plot(p2)
      dev.off()
      
      
    }
    
  }
  
}
  

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))
