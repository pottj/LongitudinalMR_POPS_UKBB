#' ---
#' title: "Get Forest Plots"
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
#' I want to repeat this for all my sensitivity analyses. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../../SourceFile_HPC.R")

#' get to do list
myFiles = list.files(path ="../results/",pattern = "03_MVMR")
myFiles = myFiles[c(1,6)]

dumTab = foreach(i = 1:length(myFiles))%do%{
  #i=1
  mySetting = gsub(".*_","",myFiles[i])
  mySetting = gsub(".RData","",mySetting)
  
  # get data
  load(paste0("../results/",myFiles[i]))
  #MVMR_results = MVMR_results[setting == "multivariate"]
  MVMR_results = MVMR_results[grepl("all_SNPs",threshold),]
  MVMR_results[exposure_type == "slope", beta_IVW := beta_IVW/55.654]
  MVMR_results[exposure_type == "slope", SE_IVW := SE_IVW/55.654]
  
  # men vs women (2SMR)
  MVMR2 = copy(MVMR_results)
  MVMR2 = MVMR2[outcome %in% c("Agaram_CAD_combined"),]
  MVMR2 = MVMR2[!grepl("combined",exposure),]
  MVMR2[exposure == "TC_men",ID:="2SMR - men"]
  MVMR2[exposure == "TC_women",ID:="2SMR - women"]
  setorder(MVMR2,ID)
  
  # plot 2
  # maybe change colors --> a bit brighter (black text vs full color might be tricky to read) or change text color? 
  MVMR2[,rank := 2]
  dummy = data.table(exposure_type = unique(MVMR2$exposure_type),rank=1)
  data2 = rbind(MVMR2,dummy, fill=T)
  data2[,subgroup := paste0("   ", ID, " - " ,setting)]
  data2[is.na(setting),subgroup := exposure_type]
  data2[,subgroup := gsub("2SMR - ","",subgroup)]
  setnames(data2,"subgroup", "   sample setting - MR model")
  data2[,lowerCI95 := beta_IVW-1.96*SE_IVW]
  data2[,upperCI95 := beta_IVW+1.96*SE_IVW]
  data2$` ` <- paste(rep(" ", 50), collapse = " ")
  data2$`Estimate [95% CI]` <- ifelse(is.na(data2$SE_IVW), "",
                                      sprintf("%.2f [%.2f, %.2f]",
                                              data2$beta_IVW, data2$lowerCI95, data2$upperCI95))
  setorder(data2,exposure_type,rank)
  
  dummy = data2$setting
  dummy2 = data2$ID
  dummy[is.na(dummy)] = "white"
  dummy[grepl("women",dummy2) & grepl("univariate",dummy)] = "darksalmon"
  dummy[grepl("women",dummy2) & grepl("multivariate",dummy)] = "firebrick"
  dummy[grepl("men",dummy2) & grepl("univariate",dummy)] = "#DEEBF7"
  dummy[grepl("men",dummy2) & grepl("multivariate",dummy)] = "#5B9BD5"
  tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))
  
  myXlab = paste0("Causal effect of TC on CAD by exposure type (",mySetting,")")
  
  p2<- forest(data2[,c(16,19,20)],
              est = data2$beta_IVW,
              lower = data2$lowerCI95, 
              upper = data2$upperCI95,
              sizes = 0.5,
              ticks_digits = 2,
              ci_column = 2,
              ref_line = 0,
              title = myXlab,
              theme = tm1)
  
  plot(p2)
  
  filename = paste0("../results/_figures/04_ForestPlots_Slides/",mySetting,"_sexComparison_2SMR.png")
  png(filename = filename,width = 1800, height = 1000, res=200)
  plot(p2)
  dev.off()
  
}

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))
