#' ---
#' title: "Get Sensitivity Forest Plots"
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
#' **Forest Plots with sensitivity data**
#' 
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
load("../results/04_MVMR_01_MAIN.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
matchingTable = data.table(old = unique(c(MVMR_results$exposure,MVMR_results$outcome)),
                           new = c("EFW_R","EFW_L","EFW_Z","EFW_C", "BW_R","BW_C","BW_Z","eCS"))

matched = match(MVMR_results$exposure,matchingTable$old)
MVMR_results[,exposure := matchingTable[matched,new]]
matched = match(MVMR_results$outcome,matchingTable$old)
MVMR_results[,outcome := matchingTable[matched,new]]
MVMR_1Sample = copy(MVMR_results)
MVMR_1Sample[outcome == "BW_R",ID:="1SMR - BW (raw)"]
MVMR_1Sample[outcome == "BW_C",ID:="1SMR - BW (percentiles)"]
MVMR_1Sample[outcome == "BW_Z",ID:="1SMR - BW (Z-scores)"]
MVMR_1Sample[outcome == "eCS",ID:="1SMR - emergency CS"]
MVMR_1Sample = MVMR_1Sample[outcome == "eCS",]
MVMR_1Sample = MVMR_1Sample[exposure == "EFW_L",]
MVMR_1Sample = MVMR_1Sample[grepl("SNPs",threshold),]

load("../results/07_MVMR_2SampleMR_CS.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
matched = match(MVMR_results$exposure,matchingTable$old)
MVMR_results[,exposure := matchingTable[matched,new]]
MVMR_2sample_BW = copy(MVMR_results)
MVMR_2sample_BW[outcome=="CS_Sakaue",ID:="2SMR - CS (Sakaue et al.)"]
MVMR_2sample_BW[outcome=="CS_UKB_elective ",ID:="2SMR - CS (elective) UKBB"]
MVMR_2sample_BW[outcome=="CS_UKB_emergency",ID:="2SMR - CS (emergency) UKBB"]
MVMR_2sample_BW = MVMR_2sample_BW[outcome == "CS_Sakaue",]
MVMR_2sample_BW = MVMR_2sample_BW[exposure == "EFW_L",]
MVMR_2sample_BW = MVMR_2sample_BW[grepl("SNPs",threshold),]

MVMR = rbind(MVMR_1Sample,MVMR_2sample_BW)
MVMR[,.N,ID]

#' # Create Forest Plots ####
#' ***
MVMR[exposure_type == "slope", beta_IVW := beta_IVW/40.316]
MVMR[exposure_type == "slope", SE_IVW := SE_IVW/40.316]

# start plotting 
MVMR2 = copy(MVMR)
MVMR2[,ID := gsub(" - .*","",ID)]
MVMR2[,threshold := gsub("_"," ",threshold)]
MVMR2[,rank := 2]
dummy = data.table(exposure_type = unique(MVMR2$exposure_type),rank=1)
data4 = rbind(MVMR2,dummy, fill=T)
data4[,subgroup := paste0("   ", ID, " - " ,threshold)]
data4[is.na(threshold),subgroup := exposure_type]
setnames(data4,"subgroup", "   sample setting - SNP selection")

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

setorder(data4,exposure_type,rank)
# data4 = data4[c(15:21,8:14,1:7,22:28),]
#data4[threshold == "nominal",threshold := "p<0.05, all SNPs"]

dummy = data4$threshold
dummy2 = data4$ID
dummy[is.na(dummy)] = "white"
dummy[grepl("2SMR",dummy2) & grepl("all SNPs",dummy)] = "#E2F0D9"
dummy[grepl("1SMR",dummy2) & grepl("all SNPs",dummy)] = "#70AD47"
dummy[grepl("2SMR",dummy2) & grepl("nominal SNPs",dummy)] = "#DEEBF7"
dummy[grepl("1SMR",dummy2) & grepl("nominal SNPs",dummy)] = "#5B9BD5"
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab = 
  "Causal effect of EFW on CS by exposure type"

p2<- forest(data4[,c(16,19,20)],
            est = data4$beta_IVW,
            lower = data4$lowerCI95, 
            upper = data4$upperCI95,
            sizes = 0.5,
            #ticks_at = myTicks,ticks_digits = 2,
            ticks_at = c(-6,-4,-2,0,2,4,6,8,10,12),ticks_digits = 0,
            ci_column = 2,
            ref_line = 0,
            #xlim = c(min_data4, max_data4),
            xlim = c(-6, 12),
            title = myXlab,
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/07_2SampleMR_04_ForestPlots_vs_1SMR/EFWL_BWR_SlidesSpecial_CS.png")
png(filename = filename,width = 1800, height = 1000, res=200)
plot(p2)
dev.off()

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))
