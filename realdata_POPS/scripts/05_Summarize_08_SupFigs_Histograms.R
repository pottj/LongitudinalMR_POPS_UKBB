#' ---
#' title: "Get Supplemental Figures for real data (histogram)"
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
#' **Supplemental Figures: Histogram of observations**
#' 
#' - MAIN (all samples)
#' - SUBSET (EUR)
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../../SourceFile.R")

#' # Main data set ####
#' ***
myGD_files = list.files(path = POPS_phenotypes,pattern = "01_Prep_04")
loaded1 = load(paste0(POPS_phenotypes,myGD_files[1]))
loaded1

plotData = myTab_X[!is.na(efwcomb),.N,by=POPSID]

ggp1 = ggplot(plotData, aes(x=N)) +
  geom_histogram( binwidth=0.5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  labs(x="number of observations per individual") +
  xlim(0,4)
ggp1
ggsave(ggp1,file = "../results/_figures/SupFigs/Histogram_observations_main.png")

#' # Sensitivity check data set ####
#' ***
plotData = myTab_X[!is.na(efwcomb) & ancestry=="GBR",.N,by=POPSID]

ggp2 = ggplot(plotData[N==3,], aes(x=N)) +
  geom_histogram( binwidth=0.5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  labs(x="number of observations per individual") +
  xlim(0,4)
ggp2
ggsave(ggp2,file = "../results/_figures/SupFigs/Histogram_observations_sensitivity.png")

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))
