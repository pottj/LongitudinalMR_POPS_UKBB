#' ---
#' title: "Get Forest plots by exposure types"
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
#' **Forest Plots of the main analysis**
#' 
#' MVMR results for EFW on BW and eCS in the main analysis. 
#' 
#' I want for each exposure type (mean, slope, var) and outcome (BW_R and eCS) one Forest Plot over all exposures (EFW_R, EFW_L, EFW_Z, EFW_C).
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
load("../results/07_MVMR_2SampleMR_BW.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
matchingTable = data.table(old = unique(c(MVMR_results$exposure)),
                           new = c("EFW_R","EFW_L","EFW_Z","EFW_C"))
matched = match(MVMR_results$exposure,matchingTable$old)
MVMR_results[,exposure := matchingTable[matched,new]]
MVMR_BW = copy(MVMR_results)

load("../results/07_MVMR_2SampleMR_CS.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
matchingTable = data.table(old = unique(c(MVMR_results$exposure)),
                           new = c("EFW_R","EFW_L","EFW_Z","EFW_C"))
matched = match(MVMR_results$exposure,matchingTable$old)
MVMR_results[,exposure := matchingTable[matched,new]]
MVMR_CS = copy(MVMR_results)

MVMR_results = rbind(MVMR_BW,MVMR_CS)

#' # Create Forest Plots ####
#' ***
#' ## BW - Mean ####
MVMR3 = copy(MVMR_results)
MVMR3 = MVMR3[exposure_type == "mean" & outcome == "BW_UKBraw",]
MVMR3[,rank := 2]

# add some more rows for the subgroup headers
dummy = data.table(exposure = unique(MVMR3$exposure),rank=1)
data4 = rbind(MVMR3,dummy, fill=T)
data4[,subgroup := paste0("   ", threshold)]
data4[is.na(threshold),subgroup := exposure]
data4[,subgroup := gsub("_", " - ",subgroup)]
setnames(data4,"subgroup", "   Instrument selection")

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
setorder(data4,exposure,rank)

dummy = data4$threshold
dummy[is.na(dummy)] = "white"
dummy[grepl("all_SNPs",dummy)] = "lightgrey"
dummy[grepl("nominal_SNPs",dummy)] = "lightsteelblue2"
dummy = gsub("top20_overlap","lightgreen",dummy)
dummy = gsub("top20_distinct","lightcoral",dummy)
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab = paste0("Mean of EFW on birth weight (UKB raw)")

p2<- forest(data4[,c(16,19,20)],
            est = data4$beta_IVW,
            lower = data4$lowerCI95, 
            upper = data4$upperCI95,
            sizes = 0.5,
            ticks_at = c(-1,0,2.5,5,7.5,10),ticks_digits = 1,
            ci_column = 2,
            ref_line = 0,
            xlim = c(-1, 10),
            title = myXlab,
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/07_2SampleMR_03_ForestPlots_byExposureType/EFW_BW_Mean.png")
png(filename = filename,width = 2000, height = 1600, res=200)
plot(p2)
dev.off()

#' ## BW - Slope ####
MVMR3 = copy(MVMR_results)
MVMR3 = MVMR3[exposure_type == "slope" & outcome == "BW_UKBraw",]
MVMR3[,rank := 2]
MVMR3[,beta_IVW := beta_IVW/40]
MVMR3[,SE_IVW := SE_IVW/40]

# add some more rows for the subgroup headers
dummy = data.table(exposure = unique(MVMR3$exposure),rank=1)
data4 = rbind(MVMR3,dummy, fill=T)
data4[,subgroup := paste0("   ", threshold)]
data4[is.na(threshold),subgroup := exposure]
data4[,subgroup := gsub("_", " - ",subgroup)]
setnames(data4,"subgroup", "   Instrument selection")

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

setorder(data4,exposure,rank)

dummy = data4$threshold
dummy[is.na(dummy)] = "white"
dummy[grepl("all_SNPs",dummy)] = "lightgrey"
dummy[grepl("nominal_SNPs",dummy)] = "lightsteelblue2"
dummy = gsub("top20_overlap","lightgreen",dummy)
dummy = gsub("top20_distinct","lightcoral",dummy)
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab = paste0("Slope of EFW on birth weight (UKB raw)")

p2<- forest(data4[,c(16,19,20)],
            est = data4$beta_IVW,
            lower = data4$lowerCI95, 
            upper = data4$upperCI95,
            sizes = 0.5,
            #ticks_at = myTicks,ticks_digits = 2,
            ticks_at = c(0,1,2,3,4,5),ticks_digits = 0,
            ci_column = 2,
            ref_line = 0,
            #xlim = c(min_data4, max_data4),
            xlim = c(0, 5),
            title = myXlab,
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/07_2SampleMR_03_ForestPlots_byExposureType/EFW_BW_Slope.png")
png(filename = filename,width = 2000, height = 1600, res=200)
plot(p2)
dev.off()

#' ## BW - Variability ####
MVMR3 = copy(MVMR_results)
MVMR3 = MVMR3[exposure_type == "var" & outcome == "BW_UKBraw",]
MVMR3[,rank := 2]

# add some more rows for the subgroup headers
dummy = data.table(exposure = unique(MVMR3$exposure),rank=1)
data4 = rbind(MVMR3,dummy, fill=T)
data4[,subgroup := paste0("   ", threshold)]
data4[is.na(threshold),subgroup := exposure]
data4[,subgroup := gsub("_", " - ",subgroup)]
setnames(data4,"subgroup", "   Instrument selection")

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

setorder(data4,exposure,rank)

dummy = data4$threshold
dummy[is.na(dummy)] = "white"
dummy[grepl("all",dummy)] = "lightgrey"
dummy[grepl("nominal",dummy)] = "lightsteelblue2"
dummy = gsub("top20_overlap","lightgreen",dummy)
dummy = gsub("top20_distinct","lightcoral",dummy)
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab = paste0("Variability of EFW on birth weight (UKB raw)")

p2<- forest(data4[,c(16,19,20)],
            est = data4$beta_IVW,
            lower = data4$lowerCI95, 
            upper = data4$upperCI95,
            sizes = 0.5,
            #ticks_at = myTicks,ticks_digits = 2,
            ticks_at = c(-0.1,-0,0.1,0.2,0.3),ticks_digits = 2,
            ci_column = 2,
            ref_line = 0,
            #xlim = c(min_data4, max_data4),
            xlim = c(-0.1, 0.35),
            title = myXlab,
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/07_2SampleMR_03_ForestPlots_byExposureType/EFW_BW_Variability.png")
png(filename = filename,width = 2000, height = 1600, res=200)
plot(p2)
dev.off()

#' ## CS - Mean ####
MVMR3 = copy(MVMR_results)
MVMR3 = MVMR3[exposure_type == "mean" & outcome == "CS_Sakaue",]
MVMR3[,rank := 2]

# add some more rows for the subgroup headers
dummy = data.table(exposure = unique(MVMR3$exposure),rank=1)
data4 = rbind(MVMR3,dummy, fill=T)
data4[,subgroup := paste0("   ", threshold)]
data4[is.na(threshold),subgroup := exposure]
data4[,subgroup := gsub("_", " - ",subgroup)]
setnames(data4,"subgroup", "   Instrument selection")

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
setorder(data4,exposure,rank)

dummy = data4$threshold
dummy[is.na(dummy)] = "white"
dummy[grepl("all_SNPs",dummy)] = "lightgrey"
dummy[grepl("nominal_SNPs",dummy)] = "lightsteelblue2"
dummy = gsub("top20_overlap","lightgreen",dummy)
dummy = gsub("top20_distinct","lightcoral",dummy)
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab = paste0("Mean of EFW on Caesarean Section (Sakaue et al.)")

p2<- forest(data4[,c(16,19,20)],
            est = data4$beta_IVW,
            lower = data4$lowerCI95, 
            upper = data4$upperCI95,
            sizes = 0.5,
            #ticks_at = myTicks,ticks_digits = 2,
            ticks_at = c(-5,-2.5,0,2.5,5),ticks_digits = 1,
            ci_column = 2,
            ref_line = 0,
            #xlim = c(min_data4, max_data4),
            xlim = c(-5, 5),
            title = myXlab,
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/07_2SampleMR_03_ForestPlots_byExposureType/EFW_CS_Mean.png")
png(filename = filename,width = 2000, height = 1600, res=200)
plot(p2)
dev.off()

#' ## CS - Slope ####
MVMR3 = copy(MVMR_results)
MVMR3 = MVMR3[exposure_type == "slope" & outcome == "CS_Sakaue",]
MVMR3[,rank := 2]
MVMR3[,beta_IVW := beta_IVW/40]
MVMR3[,SE_IVW := SE_IVW/40]

# add some more rows for the subgroup headers
dummy = data.table(exposure = unique(MVMR3$exposure),rank=1)
data4 = rbind(MVMR3,dummy, fill=T)
data4[,subgroup := paste0("   ", threshold)]
data4[is.na(threshold),subgroup := exposure]
data4[,subgroup := gsub("_", " - ",subgroup)]
setnames(data4,"subgroup", "   Instrument selection")

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

setorder(data4,exposure,rank)

dummy = data4$threshold
dummy[is.na(dummy)] = "white"
dummy[grepl("all_SNPs",dummy)] = "lightgrey"
dummy[grepl("nominal_SNPs",dummy)] = "lightsteelblue2"
dummy = gsub("top20_overlap","lightgreen",dummy)
dummy = gsub("top20_distinct","lightcoral",dummy)
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab = paste0("Slope of EFW on Caesarean Section (Sakaue et al.)")

p2<- forest(data4[,c(16,19,20)],
            est = data4$beta_IVW,
            lower = data4$lowerCI95, 
            upper = data4$upperCI95,
            sizes = 0.5,
            #ticks_at = myTicks,ticks_digits = 2,
            ticks_at = c(-5,-2.5,0,2.5,5),ticks_digits = 1,
            ci_column = 2,
            ref_line = 0,
            #xlim = c(min_data4, max_data4),
            xlim = c(-5, 5),
            title = myXlab,
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/07_2SampleMR_03_ForestPlots_byExposureType/EFW_CS_Slope.png")
png(filename = filename,width = 2000, height = 1600, res=200)
plot(p2)
dev.off()

#' ## CS - Variability ####
MVMR3 = copy(MVMR_results)
MVMR3 = MVMR3[exposure_type == "var" & outcome == "CS_Sakaue",]
MVMR3[,rank := 2]

# add some more rows for the subgroup headers
dummy = data.table(exposure = unique(MVMR3$exposure),rank=1)
data4 = rbind(MVMR3,dummy, fill=T)
data4[,subgroup := paste0("   ", threshold)]
data4[is.na(threshold),subgroup := exposure]
data4[,subgroup := gsub("_", " - ",subgroup)]
setnames(data4,"subgroup", "   Instrument selection")

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

setorder(data4,exposure,rank)

dummy = data4$threshold
dummy[is.na(dummy)] = "white"
dummy[grepl("all",dummy)] = "lightgrey"
dummy[grepl("nominal",dummy)] = "lightsteelblue2"
dummy = gsub("top20_overlap","lightgreen",dummy)
dummy = gsub("top20_distinct","lightcoral",dummy)
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab = paste0("Variability of EFW on Caesarean Section (Sakaue et al.)")

p2<- forest(data4[,c(16,19,20)],
            est = data4$beta_IVW,
            lower = data4$lowerCI95, 
            upper = data4$upperCI95,
            sizes = 0.5,
            #ticks_at = myTicks,ticks_digits = 2,
            ticks_at = c(-0.5,-0.25,0,0.25,0.5),ticks_digits = 2,
            ci_column = 2,
            ref_line = 0,
            #xlim = c(min_data4, max_data4),
            xlim = c(-0.5, 0.5),
            title = myXlab,
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/07_2SampleMR_03_ForestPlots_byExposureType/EFW_CS_Variability.png")
png(filename = filename,width = 2000, height = 1600, res=200)
plot(p2)
dev.off()

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))
