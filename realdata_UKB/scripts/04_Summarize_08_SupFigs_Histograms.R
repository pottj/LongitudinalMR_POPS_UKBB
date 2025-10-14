#' ---
#' title: "Get Supplemental Figures for real data (histogram)"
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
#' **Supplemental Figures: Histogram of observations**
#' 
#' - MAIN (all samples)
#' - SUBSET (no statins, ...)
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../../SourceFile.R")

#' # Main data set ####
#' ***
load(paste0(UKB_phenotypes_filtered,"/01_Prep_01_BL_FU_GP_merged_filtered.RData"))
plotData = copy(myTab_long)[,.N,by=ID]

ggp1 = ggplot(plotData, aes(x=N)) +
  geom_histogram( binwidth=0.5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  labs(x="number of observations per individual") 
ggp1
ggsave(ggp1,file = "../results/_figures/SupFigs/Histogram_observations_main.png")

#' # Sensitivity check data set ####
#' ***
plotData = copy(myTab_long[sens1==T,])[,.N,by=ID]

ggp2 = ggplot(plotData, aes(x=N)) +
  geom_histogram( binwidth=0.5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  labs(x="number of observations per individual") 
ggp2
ggsave(ggp2,file = "../results/_figures/SupFigs/Histogram_observations_sensitivity.png")

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))
