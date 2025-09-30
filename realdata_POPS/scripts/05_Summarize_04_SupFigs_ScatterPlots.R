#' ---
#' title: "Get Supplemental Figures for real data (scatter plots)"
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
#' **Supplemental Figures: Scatter Plots**
#' 
#' - MAIN
#' - 2-sample MVMR
#' - "all" SNPs
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../../SourceFile_HPC.R")
source("../../helperfunctions/MVMR_jp_POPS_plots.R")

#' # Load data ####
#' ***
load("../temp/04_MVMRInput_MAIN.RData")
load("../results/01_Prep_02_LD_filtered_EGG.RData")

#' # Run MVMR ####
#' ***
myAssocs_Y = myAssocs_Y[phenotype == "UKB_BW",]
myAssocs_Y[,phenotype := "BW (data from UKB)",]
myAssocs_X_long = myAssocs_X_long[phenotype == "logefwcomb",]
myAssocs_X_long[,phenotype := "EFW (data from POPS)",]
myAssocs_X_long[,dumID := "plotting"]
myAssocs_X_long[type == "var",type := "variability"]

MVMR0 = MVMR_jp_POPS_plots(data_exposure = myAssocs_X_long,
                     data_outcome = myAssocs_Y,
                     exposure_name = "EFW (data from POPS)", 
                     outcome_name = "BW (data from UKB)",
                     flag = "plotting",
                     GX_pval_treshold = 1,
                     getPlot = T,
                     corTab = LDTab2,
                     corTab_threshold = 0.1,sampleSize_GX = 2996,random = F,getCondF = T,getUni = T,
                     filename1 = "../results/_figures/SupFigs/ScatterPlot_MAIN_all_2Sample.png")

MVMR1 = MVMR_jp_POPS_plots(data_exposure = myAssocs_X_long,
                           data_outcome = myAssocs_Y,
                           exposure_name = "EFW (data from POPS)", 
                           outcome_name = "BW (data from UKB)",
                           flag = "plotting",
                           GX_pval_treshold = 5e-2,
                           getPlot = T,
                           corTab = LDTab2,
                           corTab_threshold = 0.1,sampleSize_GX = 73778,random = F,getCondF = T,getUni = T,
                           filename1 = "../results/_figures/SupFigs/ScatterPlot_MAIN_nominal_2Sample.png")

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))


