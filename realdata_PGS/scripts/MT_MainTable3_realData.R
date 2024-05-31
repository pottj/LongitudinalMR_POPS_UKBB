#' ---
#' title: "Get Main Table 3"
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
#' **Main Table for real data**
#' 
#' MVMR results for EFW on BW and eCS in the main analysis. 
#' 
#' I want one table for no p-value threshold, and adding the condtitional F-Statistics from the random setting.
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../../SourceFile_HPC.R")

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("2024-","24-",tag)
tag = gsub("-","",tag)

#' # Get data ####
#' ***
load("../results/04_MVMR_01_MAIN_240414.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
matchingTable = data.table(old = unique(c(MVMR_results$exposure,MVMR_results$outcome)),
                           new = c("EFW_R","EFW_L","EFW_Z","EFW_C",
                                   "BW_R","BW_C","BW_Z","eCS"))

matched = match(MVMR_results$exposure,matchingTable$old)
MVMR_results[,exposure := matchingTable[matched,new]]
matched = match(MVMR_results$outcome,matchingTable$old)
MVMR_results[,outcome := matchingTable[matched,new]]

MVMR1 = MVMR_results[grepl("no p-value filter",threshold)]

#' # Create Table ####
#' ***
names(MVMR1)
MVMR12 = dcast(MVMR1, exposure + outcome + threshold ~ exposure_type, value.var=c("beta_IVW","SE_IVW","condFstat"))

MVMR12 = MVMR12[,c(1,2, 4,7,10, 5,8,11, 6,9,12, 3)]
names(MVMR12) = gsub("IVW_","",names(MVMR12))
MVMR2 = copy(MVMR12)
MVMR2 = MVMR2[is.na(condFstat_mean),]
MVMR2[,condFstat_mean := MVMR12[!is.na(condFstat_mean),condFstat_mean]]
MVMR2[,condFstat_slope := MVMR12[!is.na(condFstat_slope),condFstat_slope]]
MVMR2[,condFstat_var := MVMR12[!is.na(condFstat_var),condFstat_var]]

#' create some order
setorder(MVMR2,exposure,outcome)
MVMR2[,threshold:=NULL]

MVMR2[,beta_mean := round(beta_mean,2)]
MVMR2[,beta_slope := round(beta_slope,2)]
MVMR2[,beta_var := round(beta_var,2)]

MVMR2[,SE_mean := round(SE_mean,2)]
MVMR2[,SE_slope := round(SE_slope,2)]
MVMR2[,SE_var := round(SE_var,2)]

MVMR2[,condFstat_mean := round(condFstat_mean,2)]
MVMR2[,condFstat_slope := round(condFstat_slope,2)]
MVMR2[,condFstat_var := round(condFstat_var,2)]

MVMR2

#' # Save ####
#' ***
print(xtable(MVMR2, type = "latex",digits = 3),file = "../results/_tables/MainTable3.tex")

excel_fn = paste0("../results/_tables/MainTable3.xlsx")
WriteXLS(MVMR2, 
         ExcelFileName=excel_fn, 
         SheetNames="Table3", 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

save(MVMR2, file = paste0("../results/_tables/MainTable3.RData"))

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))
