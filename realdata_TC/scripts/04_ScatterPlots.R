#' ---
#' title: "Get Scatter Plots"
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
source("../../helperfunctions/MVMR_jp_POPS_plots.R")

#' # MAIN ####
#' ***
load(file = paste0("../temp/03_MVMRInput_MAIN.RData"))
load("../results/01_Prep_03_LD_filtered.RData")
LDTab = copy(LDTab2)
LDTab[,SNP2 := as.character(SNP2)]

myExposures = unique(myAssocs_X_long$model)
mySampleSize = c(73778, 35726, 38052)
myOutcomes = unique(myAssocs_Y$phenotype)
setnames(myAssocs_Y,"SNP","markername")
setnames(myAssocs_Y,"rsID","SNP")

myAssocs_X_long[,phenotype := "Total Cholesterol"]
myAssocs_X_long[,dumID := "plotting"]
myAssocs_Y[phenotype == "Agaram_CAD_combined",phenotype := "CAD"]

filename = paste0("../results/_figures/04_ScatterPlots/MAIN_sexCombined.png")
png(filename = filename,width = 1800, height = 1000, res=200)
MVMR0 = MVMR_jp_POPS_plots(data_exposure = myAssocs_X_long[model == "combined",],
                           data_outcome = myAssocs_Y[phenotype == "CAD"],
                           exposure_name = "Total Cholesterol", 
                           outcome_name = "CAD",
                           flag = "combined",
                           GX_pval_treshold = 1e-6,
                           getPlot = T,
                           corTab = LDTab,
                           corTab_threshold = 0.1,
                           sampleSize_GX = mySampleSize[1],
                           random = F,getCondF = T,getUni = T)
dev.off()

filename = paste0("../results/_figures/04_ScatterPlots/MAIN_men.png")
png(filename = filename,width = 1800, height = 1000, res=200)
MVMR1 = MVMR_jp_POPS_plots(data_exposure = myAssocs_X_long[model == "men",],
                           data_outcome = myAssocs_Y[phenotype == "CAD"],
                           exposure_name = "Total Cholesterol", 
                           outcome_name = "CAD",
                           flag = "men",
                           GX_pval_treshold = 1e-6,
                           getPlot = T,
                           corTab = LDTab,
                           corTab_threshold = 0.1,
                           sampleSize_GX = mySampleSize[1],
                           random = F,getCondF = T,getUni = T)
dev.off()

filename = paste0("../results/_figures/04_ScatterPlots/MAIN_women.png")
png(filename = filename,width = 1800, height = 1000, res=200)
MVMR2 = MVMR_jp_POPS_plots(data_exposure = myAssocs_X_long[model == "women",],
                           data_outcome = myAssocs_Y[phenotype == "CAD"],
                           exposure_name = "Total Cholesterol", 
                           outcome_name = "CAD",
                           flag = "women",
                           GX_pval_treshold = 1e-6,
                           getPlot = T,
                           corTab = LDTab,
                           corTab_threshold = 0.1,
                           sampleSize_GX = mySampleSize[1],
                           random = F,getCondF = T,getUni = T)
dev.off()

MVMR_main = rbind(MVMR0,MVMR1,MVMR2)
MVMR_main[,flag := "main"]

#' # SENS ####
#' ***
load(file = paste0("../temp/03_MVMRInput_SENS_SNPset2.RData"))
load("../results/01_Prep_06_LD_filtered.RData")
LDTab = copy(LDTab2)
LDTab[,SNP2 := as.character(SNP2)]

myExposures = unique(myAssocs_X_long$model)
mySampleSize = c(73778, 35726, 38052)
myOutcomes = unique(myAssocs_Y$phenotype)
setnames(myAssocs_Y,"SNP","markername")
setnames(myAssocs_Y,"rsID","SNP")

myAssocs_X_long[,phenotype := "Total Cholesterol"]
myAssocs_X_long[,dumID := "plotting"]
myAssocs_Y[phenotype == "Agaram_CAD_combined",phenotype := "CAD"]

filename = paste0("../results/_figures/04_ScatterPlots/SENS_sexCombined.png")
png(filename = filename,width = 1800, height = 1000, res=200)
MVMR0 = MVMR_jp_POPS_plots(data_exposure = myAssocs_X_long[model == "combined",],
                           data_outcome = myAssocs_Y[phenotype == "CAD"],
                           exposure_name = "Total Cholesterol", 
                           outcome_name = "CAD",
                           flag = "combined",
                           GX_pval_treshold = 1e-6,
                           getPlot = T,
                           corTab = LDTab,
                           corTab_threshold = 0.1,
                           sampleSize_GX = mySampleSize[1],
                           random = F,getCondF = T,getUni = T)
dev.off()

filename = paste0("../results/_figures/04_ScatterPlots/SENS_men.png")
png(filename = filename,width = 1800, height = 1000, res=200)
MVMR1 = MVMR_jp_POPS_plots(data_exposure = myAssocs_X_long[model == "men",],
                           data_outcome = myAssocs_Y[phenotype == "CAD"],
                           exposure_name = "Total Cholesterol", 
                           outcome_name = "CAD",
                           flag = "men",
                           GX_pval_treshold = 1e-6,
                           getPlot = T,
                           corTab = LDTab,
                           corTab_threshold = 0.1,
                           sampleSize_GX = mySampleSize[1],
                           random = F,getCondF = T,getUni = T)
dev.off()

filename = paste0("../results/_figures/04_ScatterPlots/SENS_women.png")
png(filename = filename,width = 1800, height = 1000, res=200)
MVMR2 = MVMR_jp_POPS_plots(data_exposure = myAssocs_X_long[model == "women",],
                           data_outcome = myAssocs_Y[phenotype == "CAD"],
                           exposure_name = "Total Cholesterol", 
                           outcome_name = "CAD",
                           flag = "women",
                           GX_pval_treshold = 1e-6,
                           getPlot = T,
                           corTab = LDTab,
                           corTab_threshold = 0.1,
                           sampleSize_GX = mySampleSize[1],
                           random = F,getCondF = T,getUni = T)
dev.off()

MVMR_sens = rbind(MVMR0,MVMR1,MVMR2)
MVMR_sens[,flag := "sens"]

MVMR = rbind(MVMR_main,MVMR_sens)

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

