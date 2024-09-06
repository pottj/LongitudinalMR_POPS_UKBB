#' ---
#' title: "Get Forest Plots"
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

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("2024-","24-",tag)
tag = gsub("-","",tag)

#' # TC combined on CAD ####
#' ***
load("../results/03_MVMR_01_MAIN.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
MVMR_results = MVMR_results[grepl("combined",exposure),]
MVMR_results = MVMR_results[grepl("combined",outcome),]
MVMR_results[outcome == "FinnGen_UKB_CAD_combined",ID:="1SMR (UKBB for FinnGen meta)"]
#MVMR_results[outcome == "UKB_CAD_combined",ID:="1SMR (UKBB in Neale lab)"]
MVMR_results[outcome == "Agaram_CAD_combined",ID:="2SMR (Aragam et al.)"]
MVMR_results = MVMR_results[grepl("SNPs",threshold),]

MVMR = copy(MVMR_results)
MVMR = MVMR[ID != "main"]
MVMR[,.N,ID]

MVMR[exposure_type == "slope", beta_IVW := beta_IVW/55.654]
MVMR[exposure_type == "slope", SE_IVW := SE_IVW/55.654]
setorder(MVMR,ID)

# start plotting 
MVMR2 = copy(MVMR)
#MVMR2[,ID := gsub(" - .*","",ID)]
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
#dummy[grepl("1SMR",dummy2) & grepl("FinnGen",dummy2) & grepl("all SNPs",dummy)] = "darkgreen"
dummy[grepl("1SMR",dummy2) & grepl("all SNPs",dummy)] = "#70AD47"
dummy[grepl("2SMR",dummy2) & grepl("nominal SNPs",dummy)] = "#DEEBF7"
#dummy[grepl("1SMR",dummy2) & grepl("FinnGen",dummy2) & grepl("nominal SNPs",dummy)] = "darkblue"
dummy[grepl("1SMR",dummy2) & grepl("nominal SNPs",dummy)] = "#5B9BD5"
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab = 
  "Causal effect of TC on CAD by exposure type"

p2<- forest(data4[,c(16,19,20)],
            est = data4$beta_IVW,
            lower = data4$lowerCI95, 
            upper = data4$upperCI95,
            sizes = 0.5,
            #ticks_at = myTicks,ticks_digits = 2,
            #ticks_at = c(0,1,2,3,4,5,6),ticks_digits = 0,
            ci_column = 2,
            ref_line = 0,
            #xlim = c(min_data4, max_data4),
            #xlim = c(-0.25, 6),
            title = myXlab,
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/04_ForestPlots_combined.png")
png(filename = filename,width = 1800, height = 1000, res=200)
plot(p2)
dev.off()

#' # TC on CAD, men vs women (1SMR) ####
#' ***
load("../results/03_MVMR_01_MAIN.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
MVMR_results = MVMR_results[(grepl("_men",exposure) & outcome %in% c("UKB_CAD_males")) | (grepl("_women",exposure) & outcome %in% c("UKB_CAD_females")),]
MVMR_results[outcome == "UKB_CAD_males",ID:="1SMR - men"]
MVMR_results[outcome == "UKB_CAD_females",ID:="1SMR - women"]
MVMR_results = MVMR_results[grepl("SNPs",threshold),]

MVMR = copy(MVMR_results)
MVMR[,.N,ID]

MVMR[exposure_type == "slope", beta_IVW := beta_IVW/55.654]
MVMR[exposure_type == "slope", SE_IVW := SE_IVW/55.654]
setorder(MVMR,ID)

# start plotting 
MVMR2 = copy(MVMR)
#MVMR2[,ID := gsub(" - .*","",ID)]
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
dummy[grepl(" - women",dummy2) & grepl("all SNPs",dummy)] = "darksalmon"
dummy[grepl(" - women",dummy2) & grepl("nominal SNPs",dummy)] = "firebrick"
dummy[grepl(" - men",dummy2) & grepl("all SNPs",dummy)] = "#DEEBF7"
dummy[grepl(" - men",dummy2) & grepl("nominal SNPs",dummy)] = "#5B9BD5"
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab = 
  "Causal effect of TC on CAD by exposure type (men vs women)"

p2<- forest(data4[,c(16,19,20)],
            est = data4$beta_IVW,
            lower = data4$lowerCI95, 
            upper = data4$upperCI95,
            sizes = 0.5,
            #ticks_at = myTicks,ticks_digits = 2,
            #ticks_at = c(0,1,2,3,4,5,6),ticks_digits = 0,
            ci_column = 2,
            ref_line = 0,
            #xlim = c(min_data4, max_data4),
            #xlim = c(-0.25, 6),
            title = myXlab,
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/04_ForestPlots_sexComparison_1SMR.png")
png(filename = filename,width = 1800, height = 1000, res=200)
plot(p2)
dev.off()

#' # TC in women on CAD ####
#' ***
load("../results/03_MVMR_01_MAIN.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
MVMR_results = MVMR_results[outcome %in% c("Agaram_CAD_combined"),]
MVMR_results = MVMR_results[!grepl("combined",exposure),]
MVMR_results[exposure == "TC_men",ID:="2SMR - men"]
MVMR_results[exposure == "TC_women",ID:="2SMR - women"]
MVMR_results = MVMR_results[grepl("SNPs",threshold),]

MVMR = copy(MVMR_results)
MVMR[,.N,ID]

MVMR[exposure_type == "slope", beta_IVW := beta_IVW/55.654]
MVMR[exposure_type == "slope", SE_IVW := SE_IVW/55.654]
setorder(MVMR,ID)

# start plotting 
MVMR2 = copy(MVMR)
#MVMR2[,ID := gsub(" - .*","",ID)]
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
dummy[grepl(" - women",dummy2) & grepl("all SNPs",dummy)] = "darksalmon"
dummy[grepl(" - women",dummy2) & grepl("nominal SNPs",dummy)] = "firebrick"
dummy[grepl(" - men",dummy2) & grepl("all SNPs",dummy)] = "#DEEBF7"
dummy[grepl(" - men",dummy2) & grepl("nominal SNPs",dummy)] = "#5B9BD5"
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab = 
  "Causal effect of TC on CAD by exposure type (men vs women)"

p2<- forest(data4[,c(16,19,20)],
            est = data4$beta_IVW,
            lower = data4$lowerCI95, 
            upper = data4$upperCI95,
            sizes = 0.5,
            #ticks_at = myTicks,ticks_digits = 2,
            #ticks_at = c(0,1,2,3,4,5,6),ticks_digits = 0,
            ci_column = 2,
            ref_line = 0,
            #xlim = c(min_data4, max_data4),
            #xlim = c(-0.25, 6),
            title = myXlab,
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/04_ForestPlots_sexComparison_2SMR.png")
png(filename = filename,width = 1800, height = 1000, res=200)
plot(p2)
dev.off()

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))
