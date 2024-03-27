#' ---
#' title: "Get Forest Plots"
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
#' **Forest Plots** for all exposures on BW (including the sensitivity runs LMM, and GBR3), seperate for mean, slope and var
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../../SourceFile_HPC.R")

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("2024-","24-",tag)
tag = gsub("-","",tag)

#' # Load data ####
#' ***
myMVMR_files = list.files(path = "../results/",pattern = "04_MVMR")
myMVMR_files = myMVMR_files[grepl("240325",myMVMR_files)]
myMVMR_files = myMVMR_files[c(1,2,5,7,8)]

dumTab1 = foreach(i=1:length(myMVMR_files))%do%{
  #i=1
  loaded1 = load(paste0("../results/",myMVMR_files[i]))
  tab = get(loaded1)
  tab
}
myTab = rbindlist(dumTab1)
myTab[,table(exposure,threshold,ID)]

#' # Filter data ####
#' ***
MVMR1 = copy(myTab)
MVMR1 = MVMR1[setting == "multivariate",]
MVMR1 = MVMR1[outcome %in% c("pn_bw"),]
MVMR1[outcome == "pn_bw", outcome := "BW"]

MVMR1[ID == "sens_GBR3", ID := "sensitivity - GBR3"]
MVMR1[ID == "sens_LMM", ID := "sensitivity - LMM"]
MVMR1[ID == "sens_it200", ID := "sensitivity - 200 iterations"]
MVMR1[ID == "sens_SigmaTimeIndep", ID := "sensitivity - no time effect on sigma"]

#' # Create Forest Plots ###
#' ***
myExposures = unique(MVMR1$exposure)

dumTab1 = foreach(i=1:length(myExposures))%do%{
  #i=4
  MVMR2 = copy(MVMR1)
  MVMR2 = MVMR2[exposure == myExposures[i],]
  
  myExposureTypes = unique(MVMR2$exposure_type)
  
  dumTab2 = foreach(j=1:length(myExposureTypes))%do%{
    #j=2
    MVMR3 = copy(MVMR2)
    MVMR3 = MVMR3[exposure_type == myExposureTypes[j],]
    MVMR3[,rank := 2]
    
    # add some more rows for the subgroup headers
    dummy = data.table(ID = unique(MVMR3$ID),rank=1)
    data4 = rbind(MVMR3,dummy, fill=T)
    data4[,subgroup := paste0("   ", threshold)]
    data4[is.na(threshold),subgroup := ID]
    data4[,subgroup := gsub("_", " - ",subgroup)]
    setnames(data4,"subgroup", "GX p-value threshold")
    
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
    
    dummy = data4$threshold
    dummy[is.na(dummy)] = "white"
    dummy = gsub("nominal","lightgrey",dummy)
    dummy = gsub("suggestive","lightsteelblue2",dummy)
    dummy = gsub("top20_overlap","lightgreen",dummy)
    dummy = gsub("top20_distinct","lightcoral",dummy)
    tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))
    
    myXlab = paste0(myExposures[i]," (",myExposureTypes[j],") on birthweight")
    
    p2<- forest(data4[,c(15,18,19)],
                est = data4$beta_IVW,
                lower = data4$lowerCI95, 
                upper = data4$upperCI95,
                sizes = 0.5,
                ticks_at = myTicks,ticks_digits = 2,
                ci_column = 2,
                ref_line = 0,
                xlim = c(min_data4, max_data4),
                title = myXlab,
                theme = tm1)
    
    plot(p2)
    
    filename = paste0("../results/_figures/05_1_ForestPlots_",myExposures[i],"_",myExposureTypes[j],"_BW.png")
    png(filename = filename,width = 2000, height = 1500, res=200)
    plot(p2)
    dev.off()
    
    filename
  }
}


#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))
