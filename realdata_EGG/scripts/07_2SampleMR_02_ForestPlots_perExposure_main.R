#' ---
#' title: "Get Forest Plots of 2Sample MR"
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
#' **Forest Plots for real data in 2sample MR setting**
#' 
#' Forest plots per exposure on each outcome (generate all, but choose only the EFW_L on BW_R)
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../../SourceFile_HPC.R")

#' # Get data ####
#' ***
load("../results/07_MVMR_2SampleMR_BW.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
matchingTable = data.table(old = unique(c(MVMR_results$exposure)),
                           new = c("EFW_R","EFW_L","EFW_Z","EFW_C"))
matched = match(MVMR_results$exposure,matchingTable$old)
MVMR_results[,exposure := matchingTable[matched,new]]
MVMR_BW = copy(MVMR_results)

load("../results/07_MVMR_2SampleMR_CS.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
matchingTable = data.table(old = unique(c(MVMR_results$exposure)),
                           new = c("EFW_R","EFW_L","EFW_Z","EFW_C"))
matched = match(MVMR_results$exposure,matchingTable$old)
MVMR_results[,exposure := matchingTable[matched,new]]
MVMR_CS = copy(MVMR_results)

MVMR_results = rbind(MVMR_BW,MVMR_CS)

#' # Create Forest Plots in Loop ####
#' ***
myOutcomes = unique(MVMR_results$outcome)
myExposures = unique(MVMR_results$exposure)

MVMR_results[exposure_type == "slope", beta_IVW := beta_IVW/40.316]
MVMR_results[exposure_type == "slope", SE_IVW := SE_IVW/40.316]

dumTab0 = foreach(k = 1:length(myOutcomes))%do%{
  #k=1
  myOutcome = myOutcomes[k]
  message("Working on outcome ",myOutcome)
  MVMR0 = copy(MVMR_results)
  MVMR0 = MVMR0[outcome == myOutcome,]
  
  dumTab1 = foreach(i = 1:length(myExposures))%do%{
    #i=1
    myExposure = myExposures[i]
    message("   and exposure ",myExposure)
    MVMR1 = copy(MVMR0)
    MVMR1 = MVMR1[exposure == myExposure]
    
    # start plotting 
    MVMR1[,rank := 2]
    dummy = data.table(exposure_type = unique(MVMR1$exposure_type),rank=1)
    data4 = rbind(MVMR1,dummy, fill=T)
    data4[,subgroup := paste0("   ", threshold)]
    data4[is.na(threshold),subgroup := exposure_type]
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
    
    setorder(data4,exposure_type,rank)
    # data4 = data4[c(15:21,8:14,1:7,22:28),]
    #data4[threshold == "nominal",threshold := "p<0.05, all SNPs"]
    
    dummy = data4$threshold
    dummy[is.na(dummy)] = "white"
    dummy[grepl("all",dummy)] = "lightgrey"
    dummy[grepl("nominal",dummy)] = "lightsteelblue2"
    dummy = gsub("top20_overlap","lightgreen",dummy)
    dummy = gsub("top20_distinct","lightcoral",dummy)
    tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))
    
    myXlab = paste0(myExposure," on ",myOutcome)
    
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
    filename = paste0("../results/_figures/07_2SampleMR_02_ForestPlots_byExposure/ForestPlots_",myX,"_",myY,".png")
    png(filename = filename,width = 2000, height = 1800, res=200)
    plot(p2)
    dev.off()
    
  }
  
}

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))
