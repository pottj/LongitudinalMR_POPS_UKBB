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
myFiles = myFiles[-1]

dumTab = foreach(i = 1:length(myFiles))%do%{
  #i=1
  mySetting = gsub(".*_","",myFiles[i])
  mySetting = gsub(".RData","",mySetting)
  
  # get data
  load(paste0("../results/",myFiles[i]))
  MVMR_results = MVMR_results[setting == "multivariate"]
  MVMR_results = MVMR_results[grepl("SNPs",threshold),]
  MVMR_results[exposure_type == "slope", beta_IVW := beta_IVW/55.654]
  MVMR_results[exposure_type == "slope", SE_IVW := SE_IVW/55.654]
  
  # sex-combined
  MVMR0 = copy(MVMR_results)
  MVMR0 = MVMR0[grepl("combined",outcome) & grepl("combined",exposure),]
  MVMR0[outcome == "FinnGen_UKB_CAD_combined",ID:="1SMR (UKBB for FinnGen meta)"]
  MVMR0[outcome == "Agaram_CAD_combined",ID:="2SMR (Aragam et al.)"]
  MVMR0 = MVMR0[!grepl("sens",ID)]
  setorder(MVMR0,ID)
  
  # men vs women (1SMR)
  MVMR1 = copy(MVMR_results)
  MVMR1 = MVMR1[(grepl("_men",exposure) & outcome %in% c("UKB_CAD_males")) | (grepl("_women",exposure) & outcome %in% c("UKB_CAD_females")),]
  MVMR1[outcome == "UKB_CAD_males",ID:="1SMR - men"]
  MVMR1[outcome == "UKB_CAD_females",ID:="1SMR - women"]
  setorder(MVMR1,ID)
  
  # men vs women (2SMR)
  MVMR2 = copy(MVMR_results)
  MVMR2 = MVMR2[outcome %in% c("Agaram_CAD_combined"),]
  MVMR2 = MVMR2[!grepl("combined",exposure),]
  MVMR2[exposure == "TC_men",ID:="2SMR - men"]
  MVMR2[exposure == "TC_women",ID:="2SMR - women"]
  setorder(MVMR2,ID)
  
  # plot 0
  MVMR0[,threshold := gsub("_"," ",threshold)]
  MVMR0[,rank := 2]
  dummy = data.table(exposure_type = unique(MVMR0$exposure_type),rank=1)
  data0 = rbind(MVMR0,dummy, fill=T)
  data0[,subgroup := paste0("   ", ID, " - " ,threshold)]
  data0[is.na(threshold),subgroup := exposure_type]
  setnames(data0,"subgroup", "   sample setting - SNP selection")
  data0[,lowerCI95 := beta_IVW-1.96*SE_IVW]
  data0[,upperCI95 := beta_IVW+1.96*SE_IVW]
  data0$` ` <- paste(rep(" ", 50), collapse = " ")
  data0$`Estimate [95% CI]` <- ifelse(is.na(data0$SE_IVW), "",
                                      sprintf("%.2f [%.2f, %.2f]",
                                              data0$beta_IVW, data0$lowerCI95, data0$upperCI95))
  setorder(data0,exposure_type,rank)
  
  dummy = data0$threshold
  dummy2 = data0$ID
  dummy[is.na(dummy)] = "white"
  dummy[grepl("2SMR",dummy2) & grepl("all SNPs",dummy)] = "#E2F0D9"
  dummy[grepl("1SMR",dummy2) & grepl("all SNPs",dummy)] = "#70AD47"
  dummy[grepl("2SMR",dummy2) & grepl("nominal SNPs",dummy)] = "#DEEBF7"
  dummy[grepl("1SMR",dummy2) & grepl("nominal SNPs",dummy)] = "#5B9BD5"
  tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))
  
  myXlab = paste0("Causal effect of TC on CAD by exposure type (",mySetting,")")
  
  p0<- forest(data0[,c(16,19,20)],
              est = data0$beta_IVW,
              lower = data0$lowerCI95, 
              upper = data0$upperCI95,
              sizes = 0.5,
              ticks_digits = 2,
              ci_column = 2,
              ref_line = 0,
              title = myXlab,
              theme = tm1)
  
  plot(p0)
  
  filename = paste0("../results/_figures/04_ForestPlots_Sensitivity/",mySetting,"_combined.png")
  png(filename = filename,width = 1800, height = 1000, res=200)
  plot(p0)
  dev.off()
  
  # plot 1
  MVMR1[,threshold := gsub("_"," ",threshold)]
  MVMR1[,rank := 2]
  dummy = data.table(exposure_type = unique(MVMR1$exposure_type),rank=1)
  data1 = rbind(MVMR1,dummy, fill=T)
  data1[,subgroup := paste0("   ", ID, " - " ,threshold)]
  data1[is.na(threshold),subgroup := exposure_type]
  setnames(data1,"subgroup", "   sample setting - SNP selection")
  data1[,lowerCI95 := beta_IVW-1.96*SE_IVW]
  data1[,upperCI95 := beta_IVW+1.96*SE_IVW]
  data1$` ` <- paste(rep(" ", 50), collapse = " ")
  data1$`Estimate [95% CI]` <- ifelse(is.na(data1$SE_IVW), "",
                                      sprintf("%.2f [%.2f, %.2f]",
                                              data1$beta_IVW, data1$lowerCI95, data1$upperCI95))
  setorder(data1,exposure_type,rank)
  
  dummy = data1$threshold
  dummy2 = data1$ID
  dummy[is.na(dummy)] = "white"
  dummy[grepl(" - women",dummy2) & grepl("all SNPs",dummy)] = "darksalmon"
  dummy[grepl(" - women",dummy2) & grepl("nominal SNPs",dummy)] = "firebrick"
  dummy[grepl(" - men",dummy2) & grepl("all SNPs",dummy)] = "#DEEBF7"
  dummy[grepl(" - men",dummy2) & grepl("nominal SNPs",dummy)] = "#5B9BD5"
  tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))
  
  myXlab = paste0("Causal effect of TC on CAD by exposure type (men vs women; ",mySetting,")")
  
  p1<- forest(data1[,c(16,19,20)],
              est = data1$beta_IVW,
              lower = data1$lowerCI95, 
              upper = data1$upperCI95,
              sizes = 0.5,
              ticks_digits = 2,
              ci_column = 2,
              ref_line = 0,
              title = myXlab,
              theme = tm1)
  
  plot(p1)
  
  filename = paste0("../results/_figures/04_ForestPlots_Sensitivity/",mySetting,"_sexComparison_1SMR.png")
  png(filename = filename,width = 1800, height = 1000, res=200)
  plot(p1)
  dev.off()
  
  # plot 2
  MVMR2[,threshold := gsub("_"," ",threshold)]
  MVMR2[,rank := 2]
  dummy = data.table(exposure_type = unique(MVMR2$exposure_type),rank=1)
  data2 = rbind(MVMR2,dummy, fill=T)
  data2[,subgroup := paste0("   ", ID, " - " ,threshold)]
  data2[is.na(threshold),subgroup := exposure_type]
  setnames(data2,"subgroup", "   sample setting - SNP selection")
  data2[,lowerCI95 := beta_IVW-1.96*SE_IVW]
  data2[,upperCI95 := beta_IVW+1.96*SE_IVW]
  data2$` ` <- paste(rep(" ", 50), collapse = " ")
  data2$`Estimate [95% CI]` <- ifelse(is.na(data2$SE_IVW), "",
                                      sprintf("%.2f [%.2f, %.2f]",
                                              data2$beta_IVW, data2$lowerCI95, data2$upperCI95))
  setorder(data2,exposure_type,rank)
  
  dummy = data2$threshold
  dummy2 = data2$ID
  dummy[is.na(dummy)] = "white"
  dummy[grepl(" - women",dummy2) & grepl("all SNPs",dummy)] = "darksalmon"
  dummy[grepl(" - women",dummy2) & grepl("nominal SNPs",dummy)] = "firebrick"
  dummy[grepl(" - men",dummy2) & grepl("all SNPs",dummy)] = "#DEEBF7"
  dummy[grepl(" - men",dummy2) & grepl("nominal SNPs",dummy)] = "#5B9BD5"
  tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))
  
  myXlab = paste0("Causal effect of TC on CAD by exposure type (men vs women (2SMR); ",mySetting,")")
  
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
  
  filename = paste0("../results/_figures/04_ForestPlots_Sensitivity/",mySetting,"_sexComparison_2SMR.png")
  png(filename = filename,width = 1800, height = 1000, res=200)
  plot(p2)
  dev.off()
  
}

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))
