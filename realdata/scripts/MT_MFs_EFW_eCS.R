#' ---
#' title: "Get Main Table"
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
#' MVMR results for EFW on BW in the main analysis. 
#' 
#' I want one table for no p-value threshold, and a forest plot per type (mean, slope, var)
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
load("../temp/04_MVMR_01_MAIN_classic_240408.RData")
MVMR_results = MVMR_results1[grepl("efw",exposure)]
MVMR_results = MVMR_results[outcome == "pn_emcsall"]
MVMR_results = MVMR_results[setting == "multivariate"]
matchingTable = data.table(old = unique(MVMR_results$exposure),
                           new = c("EFW_R","EFW_L","EFW_Z","EFW_C"))
matched = match(MVMR_results$exposure,matchingTable$old)
MVMR_results[,exposure := matchingTable[matched,new]]

MVMR1 = MVMR_results[grepl("no p-value filter",threshold)]

#' # Create Table ####
#' ***
names(MVMR1)
MVMR12 = dcast(MVMR1, exposure + NR_SNPs_total + threshold ~ exposure_type, value.var=c("beta_IVW","SE_IVW"))

MVMR12 = MVMR12[,c(1,2,4,7,5,8,6,9)]
names(MVMR12) = gsub("IVW_","",names(MVMR12))
names(MVMR12)[2] = "NR_SNPs"
MVMR12 = MVMR12[c(6,5,4,3,2,1,8,7),]

MVMR12[,beta_mean := round(beta_mean,3)]
MVMR12[,beta_slope := round(beta_slope,3)]
MVMR12[,beta_var := round(beta_var,3)]

MVMR12[,SE_mean := signif(SE_mean,3)]
MVMR12[,SE_slope := signif(SE_slope,3)]
MVMR12[,SE_var := signif(SE_var,3)]
MVMR12

print(xtable(MVMR12, type = "latex",digits = 3),file = "../results/_tables/RealData_mainTable3.tex")

#' # Create Forest Plots ####
#' ***
#' ## Mean ####
MVMR3 = copy(MVMR_results)
MVMR3 = MVMR3[exposure_type == "mean",]
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

filename = paste0("../results/_figures/MF_ForestPlots_EFW_eCS_Mean.png")
png(filename = filename,width = 2000, height = 1600, res=200)
plot(p2)
dev.off()

#' ## Slope ####
MVMR3 = copy(MVMR_results)
MVMR3 = MVMR3[exposure_type == "slope",]
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

filename = paste0("../results/_figures/MF_ForestPlots_EFW_eCS_Slope.png")
png(filename = filename,width = 2000, height = 1600, res=200)
plot(p2)
dev.off()

#' ## Variability ####
MVMR3 = copy(MVMR_results)
MVMR3 = MVMR3[exposure_type == "var",]
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

filename = paste0("../results/_figures/MF_ForestPlots_EFW_eCS_Variability.png")
png(filename = filename,width = 2000, height = 1600, res=200)
plot(p2)
dev.off()

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))
