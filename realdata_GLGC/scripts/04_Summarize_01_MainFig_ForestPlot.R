#' ---
#' title: "Get Forest Plots"
#' subtitle: "Longitudinal MVMR in UKBB"
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
#' For slides in BSU workshop talk
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../../SourceFile_HPC.R")

#' # Loop ####
#' ***
#' I only want the 2-sample MVMR
#' 
myFiles = list.files(path = "../results/", pattern = "03_MVMR")
mySettings = gsub(".RData","",myFiles)
mySettings = gsub("03_MVMR_0._","",mySettings)

for(i in 1:length(myFiles)){
  #i=1
  load(paste0("../results/",myFiles[i]))
  MVMR_results = MVMR_results[outcome == "Aragam"]
  MVMR = copy(MVMR_results)
  MVMR[,ID := setting]
  MVMR[,.N,ID]
  MVMR[exposure_type == "slope", beta_IVW := beta_IVW/55.654]
  MVMR[exposure_type == "slope", SE_IVW := SE_IVW/55.654]
  
  # start plotting 
  MVMR2 = copy(MVMR)
  MVMR2 = MVMR2[threshold == "nominal_SNPs",]
  MVMR2[,ID := gsub(" - .*","",ID)]
  MVMR2[,threshold := gsub("_"," ",threshold)]
  MVMR2[,rank := 2]
  dummy = data.table(exposure_type = unique(MVMR2$exposure_type),rank=1)
  data4 = rbind(MVMR2,dummy, fill=T)
  data4[,subgroup := paste0("   ", ID)]
  data4[is.na(threshold),subgroup := exposure_type]
  data4[subgroup=="   multivariate",subgroup:= "   MVMR"]
  data4[subgroup=="   univariate",subgroup:= "   MR"]
  #setnames(data4,"subgroup", "   model")
  
  data4[,lowerCI95 := beta_IVW-1.96*SE_IVW]
  data4[,upperCI95 := beta_IVW+1.96*SE_IVW]
  data4$` ` <- paste(rep(" ", 50), collapse = " ")
  data4$`Estimate \n[95% CI]` <- ifelse(is.na(data4$SE_IVW), "",
                                      sprintf("%.2f [%.2f, %.2f]",
                                              data4$beta_IVW, data4$lowerCI95, data4$upperCI95))
  
  setorder(data4,exposure_type,rank)
  data4[exposure_type=="var" & is.na(setting),subgroup:= "variability"]
  dummy = data4$ID
  dummy2 = data4$exposure_type
  dummy[is.na(dummy)] = "white"
  dummy[grepl("uni",dummy) & grepl("mean",dummy2)] = "#FBE3D6"
  dummy[grepl("multi",dummy) & grepl("mean",dummy2)] = "#F2AA84"
  dummy[grepl("uni",dummy) & grepl("slope",dummy2)] = "#C2F1C8"
  dummy[grepl("multi",dummy) & grepl("slope",dummy2)] = "#47D45A"
  dummy[grepl("uni",dummy) & grepl("var",dummy2)] = "#CAEEFB"
  dummy[grepl("multi",dummy) & grepl("var",dummy2)] = "#61CBF4"
  tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))
  
  myXlab =  "causal estimand of TC on CAD"
  
  data4[,condF := round(condFstat,1)]
  data4[,condF := as.character(condF)]
  data4[is.na(setting),condF := ""]
  setnames(data4,"condF","(cond) \nF-stat")
  setnames(data4,"subgroup", "exposure type \n   MR approach")
  
  p2<- forest(data4[,c(16,19,20,21)],
              est = data4$beta_IVW,
              lower = data4$lowerCI95, 
              upper = data4$upperCI95,
              sizes = 0.5,
              ci_column = 2,
              ref_line = 0,
              #title = myXlab,
              xlab = myXlab,
              theme = tm1)
  
  plot(p2)
  
  filename = paste0("../results/_figures/MainFigures/TC_CAD_",mySettings[i],".png")
  png(filename = filename,width = 1650, height = 650, res=200)
  plot(p2)
  dev.off()
  
}

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))
