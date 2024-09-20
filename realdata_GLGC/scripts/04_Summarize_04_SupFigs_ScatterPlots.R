#' ---
#' title: "Get Supplemental Figures for real data (scatter plots)"
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
load("../temp/03_MVMRInput_MAIN.RData")
load("../results/01_Prep_04_LD.RData")

#' # Run MVMR ####
#' ***
names(myAssocs_Y_long) = c("SNP", "phenotype","sampleSize","beta_mean","SE_mean","tval_mean","pval_mean" )
myAssocs_Y_long = myAssocs_Y_long[phenotype == "Aragam",]
myAssocs_Y_long[,phenotype := "CAD (data from Aragam et al.)",]
myAssocs_X_long[,phenotype := "TC (data from UKB)",]
myAssocs_X_long[,dumID := "plotting"]
myAssocs_X_long[type == "var",type := "variability"]

MVMR0 = MVMR_jp_POPS_plots(data_exposure = myAssocs_X_long,
                     data_outcome = myAssocs_Y_long,
                     exposure_name = "TC (data from UKB)", 
                     outcome_name = "CAD (data from Aragam et al.)",
                     flag = "plotting",
                     GX_pval_treshold = 1,
                     getPlot = T,
                     corTab = LDTab,
                     corTab_threshold = 0.1,sampleSize_GX = 73778,random = F,getCondF = T,getUni = T,
                     filename1 = "../results/_figures/SupFigs/ScatterPlot_MAIN_all_2Sample.png")

MVMR1 = MVMR_jp_POPS_plots(data_exposure = myAssocs_X_long,
                           data_outcome = myAssocs_Y_long,
                           exposure_name = "TC (data from UKB)", 
                           outcome_name = "CAD (data from Aragam et al.)",
                           flag = "plotting",
                           GX_pval_treshold = 5e-8,
                           getPlot = T,
                           corTab = LDTab,
                           corTab_threshold = 0.1,sampleSize_GX = 73778,random = F,getCondF = T,getUni = T,
                           filename1 = "../results/_figures/SupFigs/ScatterPlot_MAIN_nominal_2Sample.png")

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))


