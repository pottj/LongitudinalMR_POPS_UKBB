#' ---
#' title: "Get Sup Fig Type 1: Forest plots by exposure types"
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
#' **Supplemental Figures for real data**
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
load("../results/04_MVMR_01_MAIN_240414.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
matchingTable = data.table(old = unique(c(MVMR_results$exposure,MVMR_results$outcome)),
                           new = c("EFW_R","EFW_L","EFW_Z","EFW_C",
                                   "BW_R","BW_C","BW_Z","eCS"))

matched = match(MVMR_results$exposure,matchingTable$old)
MVMR_results[,exposure := matchingTable[matched,new]]
matched = match(MVMR_results$outcome,matchingTable$old)
MVMR_results[,outcome := matchingTable[matched,new]]

#' # Create Forest Plots ####
#' ***
#' ## Mean ####
MVMR3 = copy(MVMR_results)
MVMR3 = MVMR3[exposure_type == "mean" & outcome == "BW_R",]
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
data4 = data4[!grepl("suggestive",threshold),]
data4 = data4[c(15:21,8:14,1:7,22:28),]
#data4[threshold == "nominal",threshold := "p<0.05, all SNPs"]

dummy = data4$threshold
dummy[is.na(dummy)] = "white"
dummy[grepl("no p-value filter",dummy)] = "lightgrey"
dummy[grepl("p<0.05",dummy)] = "lightsteelblue2"
dummy = gsub("top20_overlap","lightgreen",dummy)
dummy = gsub("top20_distinct","lightcoral",dummy)
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab = paste0("Mean of EFW on birth weight")

p2<- forest(data4[,c(16,19,20)],
            est = data4$beta_IVW,
            lower = data4$lowerCI95, 
            upper = data4$upperCI95,
            sizes = 0.5,
            ticks_at = c(0,1,2,3),ticks_digits = 2,
            ci_column = 2,
            ref_line = 0,
            xlim = c(-0.5, 3.5),
            title = myXlab,
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/SF_SupFig7_ForestPlots_EFW_BW_Mean.png")
png(filename = filename,width = 2000, height = 1600, res=200)
plot(p2)
dev.off()

#' ## Slope ####
MVMR3 = copy(MVMR_results)
MVMR3 = MVMR3[exposure_type == "slope" & outcome == "BW_R",]
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
data4 = data4[!grepl("suggestive",threshold),]
data4 = data4[c(15:21,8:14,1:7,22:28),]
#data4[threshold == "nominal",threshold := "p<0.05, all SNPs"]

dummy = data4$threshold
dummy[is.na(dummy)] = "white"
dummy[grepl("no p-value filter",dummy)] = "lightgrey"
dummy[grepl("p<0.05",dummy)] = "lightsteelblue2"
dummy = gsub("top20_overlap","lightgreen",dummy)
dummy = gsub("top20_distinct","lightcoral",dummy)
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab = paste0("Slope of EFW on birth weight")

p2<- forest(data4[,c(16,19,20)],
            est = data4$beta_IVW,
            lower = data4$lowerCI95, 
            upper = data4$upperCI95,
            sizes = 0.5,
            #ticks_at = myTicks,ticks_digits = 2,
            ticks_at = c(0,10,20,30,40,50,60,70,80),ticks_digits = 0,
            ci_column = 2,
            ref_line = 0,
            #xlim = c(min_data4, max_data4),
            xlim = c(-10, 90),
            title = myXlab,
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/SF_SupFig8_ForestPlots_EFW_BW_Slope.png")
png(filename = filename,width = 2000, height = 1600, res=200)
plot(p2)
dev.off()

#' ## Variability ####
MVMR3 = copy(MVMR_results)
MVMR3 = MVMR3[exposure_type == "var" & outcome == "BW_R",]
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
data4 = data4[!grepl("suggestive",threshold),]
data4 = data4[c(15:21,8:14,1:7,22:28),]
#data4[threshold == "nominal",threshold := "p<0.05, all SNPs"]

dummy = data4$threshold
dummy[is.na(dummy)] = "white"
dummy[grepl("no p-value filter",dummy)] = "lightgrey"
dummy[grepl("p<0.05",dummy)] = "lightsteelblue2"
dummy = gsub("top20_overlap","lightgreen",dummy)
dummy = gsub("top20_distinct","lightcoral",dummy)
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab = paste0("Variability of EFW on birth weight")

p2<- forest(data4[,c(16,19,20)],
            est = data4$beta_IVW,
            lower = data4$lowerCI95, 
            upper = data4$upperCI95,
            sizes = 0.5,
            #ticks_at = myTicks,ticks_digits = 2,
            ticks_at = c(-0.15,-0.10,-0.05,0,0.05,0.1,0.15),ticks_digits = 2,
            ci_column = 2,
            ref_line = 0,
            #xlim = c(min_data4, max_data4),
            xlim = c(-0.2, 0.2),
            title = myXlab,
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/SF_SupFig9_ForestPlots_EFW_BW_Variability.png")
png(filename = filename,width = 2000, height = 1600, res=200)
plot(p2)
dev.off()

#' ## Mean ####
MVMR3 = copy(MVMR_results)
MVMR3 = MVMR3[exposure_type == "mean" & outcome == "eCS",]
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
data4 = data4[!grepl("suggestive",threshold),]
data4 = data4[c(15:21,8:14,1:7,22:28),]
#data4[threshold == "nominal",threshold := "p<0.05, all SNPs"]

dummy = data4$threshold
dummy[is.na(dummy)] = "white"
dummy[grepl("no p-value filter",dummy)] = "lightgrey"
dummy[grepl("p<0.05",dummy)] = "lightsteelblue2"
dummy = gsub("top20_overlap","lightgreen",dummy)
dummy = gsub("top20_distinct","lightcoral",dummy)
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab = paste0("Mean of EFW on emergency Cesarean Section")

p2<- forest(data4[,c(16,19,20)],
            est = data4$beta_IVW,
            lower = data4$lowerCI95, 
            upper = data4$upperCI95,
            sizes = 0.5,
            #ticks_at = myTicks,ticks_digits = 2,
            ticks_at = c(-25,-20,-15,-10,-5,0,5,10),ticks_digits = 0,
            ci_column = 2,
            ref_line = 0,
            #xlim = c(min_data4, max_data4),
            xlim = c(-30, 15),
            title = myXlab,
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/SF_SupFig13_ForestPlots_EFW_eCS_Mean.png")
png(filename = filename,width = 2000, height = 1600, res=200)
plot(p2)
dev.off()

#' ## Slope ####
MVMR3 = copy(MVMR_results)
MVMR3 = MVMR3[exposure_type == "slope" & outcome == "eCS",]
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
data4 = data4[!grepl("suggestive",threshold),]
data4 = data4[c(15:21,8:14,1:7,22:28),]
#data4[threshold == "nominal",threshold := "p<0.05, all SNPs"]

dummy = data4$threshold
dummy[is.na(dummy)] = "white"
dummy[grepl("no p-value filter",dummy)] = "lightgrey"
dummy[grepl("p<0.05",dummy)] = "lightsteelblue2"
dummy = gsub("top20_overlap","lightgreen",dummy)
dummy = gsub("top20_distinct","lightcoral",dummy)
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab = paste0("Slope of EFW on emergency Cesarean Section")

p2<- forest(data4[,c(16,19,20)],
            est = data4$beta_IVW,
            lower = data4$lowerCI95, 
            upper = data4$upperCI95,
            sizes = 0.5,
            #ticks_at = myTicks,ticks_digits = 2,
            ticks_at = seq(from=-500,to=200,by=100),ticks_digits = 0,
            ci_column = 2,
            ref_line = 0,
            #xlim = c(min_data4, max_data4),
            xlim = c(-550, 250),
            title = myXlab,
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/SF_SupFig14_ForestPlots_EFW_eCS_Slope.png")
png(filename = filename,width = 2000, height = 1600, res=200)
plot(p2)
dev.off()

#' ## Variability ####
MVMR3 = copy(MVMR_results)
MVMR3 = MVMR3[exposure_type == "var" & outcome == "eCS",]
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
data4 = data4[!grepl("suggestive",threshold),]
data4 = data4[c(15:21,8:14,1:7,22:28),]
#data4[threshold == "nominal",threshold := "p<0.05, all SNPs"]

dummy = data4$threshold
dummy[is.na(dummy)] = "white"
dummy[grepl("no p-value filter",dummy)] = "lightgrey"
dummy[grepl("p<0.05",dummy)] = "lightsteelblue2"
dummy = gsub("top20_overlap","lightgreen",dummy)
dummy = gsub("top20_distinct","lightcoral",dummy)
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab = paste0("Variability of EFW on emergency Cesarean Section")

p2<- forest(data4[,c(16,19,20)],
            est = data4$beta_IVW,
            lower = data4$lowerCI95, 
            upper = data4$upperCI95,
            sizes = 0.5,
            #ticks_at = myTicks,ticks_digits = 2,
            ticks_at = c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25),ticks_digits = 2,
            ci_column = 2,
            ref_line = 0,
            #xlim = c(min_data4, max_data4),
            xlim = c(-1, 1.5),
            title = myXlab,
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/SF_SupFig15_ForestPlots_EFW_eCS_Variability.png")
png(filename = filename,width = 2000, height = 1600, res=200)
plot(p2)
dev.off()

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))
