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
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../../SourceFile_HPC.R")

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("2024-","24-",tag)
tag = gsub("-","",tag)

#' # Load and filter
#' ***
load("../results/03_MVMR_01_MAIN.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
filt0 = grepl("combined",MVMR_results$exposure) & grepl("combined",MVMR_results$outcome)
filt1 = grepl("_men",MVMR_results$exposure) & MVMR_results$outcome %in% c("Agaram_CAD_combined","FinnGen_UKB_CAD_combined","UKB_CAD_males")
filt2 = grepl("_women",MVMR_results$exposure) & MVMR_results$outcome %in% c("Agaram_CAD_combined","FinnGen_UKB_CAD_combined","UKB_CAD_females")
filt = filt0 | filt1 | filt2
MVMR_results = MVMR_results[filt]
MVMR_results = MVMR_results[grepl("SNPs",threshold),]
MVMR_results[,dumID := paste0(exposure_type,"_2SMR")]
MVMR_results[grepl("UKB",outcome),dumID := paste0(exposure_type,"_1SMR")]
MVMR_results[grepl("FinnGen",outcome),dumID := paste0(exposure_type,"_1SMR_new")]

#' # Loop 
#' ***
#' I want 6 plots, combination of exposure type and 1- vs 2-sample MR
mySetting = unique(MVMR_results$dumID)

mySetting

for(i in 1:length(mySetting)){
  #i=1
  message("Working on setting ",mySetting[i])
  MVMR = copy(MVMR_results)
  MVMR = MVMR[dumID == mySetting[i]]
  
  MVMR[,threshold := gsub("_"," ",threshold)]
  MVMR[,rank := 2]
  dummy = data.table(exposure = unique(MVMR$exposure),rank=1)
  data4 = rbind(MVMR,dummy, fill=T)
  data4[,subgroup := paste0("         ",threshold)]
  data4[is.na(threshold),subgroup := gsub(".*_","",exposure)]
  setnames(data4,"subgroup", "   instrument selection")
  
  data4[,lowerCI95 := beta_IVW-1.96*SE_IVW]
  data4[,upperCI95 := beta_IVW+1.96*SE_IVW]
  data4$` ` <- paste(rep(" ", 50), collapse = " ")
  data4$`Estimate [95% CI]` <- ifelse(is.na(data4$SE_IVW), "",
                                      sprintf("%.2f [%.2f, %.2f]",
                                              data4$beta_IVW, data4$lowerCI95, data4$upperCI95))
  
  min_data4 = round(min(data4$lowerCI95,na.rm = T) - 0.1,2)
  max_data4 = round(max(data4$upperCI95,na.rm = T) + 0.1,2)
  range = max_data4 - min_data4
  myTicks = c(seq(min_data4, max_data4, by = round(range/5,2)),max_data4)
  
  setorder(data4,exposure,rank)
  # data4 = data4[c(15:21,8:14,1:7,22:28),]
  #data4[threshold == "nominal",threshold := "p<0.05, all SNPs"]
  
  dummy = data4$threshold
  dummy2 = data4$exposure
  dummy[is.na(dummy)] = "white"
  dummy[grepl("combined",dummy2) & grepl("nominal SNPs",dummy)] = "#E2F0D9"
  dummy[grepl("combined",dummy2) & grepl("all SNPs",dummy)] = "#70AD47"
  dummy[grepl("_men",dummy2) & grepl("nominal SNPs",dummy)] = "#DEEBF7"
  dummy[grepl("_men",dummy2) & grepl("all SNPs",dummy)] = "#5B9BD5"
  dummy[grepl("_women",dummy2) & grepl("nominal SNPs",dummy)] = "darksalmon"
  dummy[grepl("_women",dummy2) & grepl("all SNPs",dummy)] = "firebrick"
  tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))
  
  myXlab = paste0("Causal effect of ",gsub("_.*","",mySetting[i])," of TC on CAD by sex (",gsub(".*_","",mySetting[i]),")")
  
  p2<- forest(data4[,c(17,20,21)],
              est = data4$beta_IVW,
              lower = data4$lowerCI95, 
              upper = data4$upperCI95,
              sizes = 0.5,
              #ticks_at = myTicks,ticks_digits = 2,
              #ticks_at = c(0,1,2,3,4,5,6),ticks_digits = 0,
              ci_column = 2,
              ref_line = 0,
              #xlim = c(min_data4, max_data4),
              #xlim = c(-0.25, 6),
              title = myXlab,
              theme = tm1)
  
  plot(p2)
  
  filename = paste0("../results/_figures/04_ForestPlots_perExposureType/",mySetting[i],".png")
  png(filename = filename,width = 1800, height = 800, res=200)
  plot(p2)
  dev.off()
  
}

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))
