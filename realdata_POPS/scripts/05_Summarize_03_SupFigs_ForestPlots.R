#' ---
#' title: "Get Supplemental Figures for real data (forest plots)"
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
#' **Supplemental Figures: Forest Plots**
#' 
#' - main & log: 2-sample (UKB) vs 1-sample (POPS)
#' - main & log: all SNPs vs SNP selection
#' - all setting: log vs Zscores & MR results
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../../SourceFile.R")

#' # Load data ####
#' ***
#' I only want the 2-sample MVMR using Aragam data and all SNPs 
#' 
myFiles = list.files(path = "../results/", pattern = "04_MVMR")
mySettings = gsub(".RData","",myFiles)
mySettings = gsub("04_MVMR_0._","",mySettings)
mySettings
mySettings2 = c("0","2","1B","1A")

dumTab2 = foreach(i = 1:length(myFiles))%do%{
  #i=1
  load(paste0("../results/",myFiles[i]))
  MVMR = copy(MVMR_results)
  MVMR[,flag := mySettings2[i]]
  MVMR
}

MVMR = rbindlist(dumTab2)

#' # Plotting ####
#' ***
#' ## 2-sample vs 1-sample ####
#' 
MVMR2 = copy(MVMR)
MVMR2 = MVMR2[threshold == "all_SNPs",]
MVMR2 = MVMR2[flag == c("0")]
MVMR2 = MVMR2[outcome != "EGG_BW",]
MVMR2 = MVMR2[exposure == "logefwcomb",]
MVMR2 = MVMR2[setting == "multivariate",]

setorder(MVMR2,flag)
MVMR2[,rank := 2]
dummy = data.table(exposure_type = unique(MVMR2$exposure_type),rank=1)
data4 = rbind(MVMR2,dummy, fill=T)
data4[,subgroup := paste0("   ", flag)]
data4[is.na(threshold),subgroup := exposure_type]

data4[,lowerCI95 := beta_IVW-1.96*SE_IVW]
data4[,upperCI95 := beta_IVW+1.96*SE_IVW]
data4$` ` <- paste(rep(" ", 50), collapse = " ")
data4$`Estimate \n[95% CI]` <- ifelse(is.na(data4$SE_IVW), "",
                                      sprintf("%.2f [%.2f, %.2f]",
                                              data4$beta_IVW, data4$lowerCI95, data4$upperCI95))

setorder(data4,exposure_type,rank)
data4[exposure_type=="var" & is.na(setting),subgroup:= "variability"]
data4[grepl("UKB",outcome), subgroup := gsub("0","main (UKB)",subgroup)]
data4[grepl("POPS",outcome), subgroup := gsub("0","sensitivity (POPS)",subgroup)]
setorder(data4,exposure_type,rank,subgroup)

dummy = c("white","#F2AA84", "#FBE3D6",
          "white","#47D45A", "#C2F1C8", 
          "white","#61CBF4", "#CAEEFB" )
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab =  "increase in BW (95% CI) per 1-SD increment EFW"

data4[,condF := round(condFstat,1)]
data4[,condF := as.character(condF)]
data4[is.na(setting),condF := ""]
setnames(data4,"condF","cond. \nF-stat")
setnames(data4,"subgroup", "Exposure type \n   Outcome data")

p2<- forest(data4[,c(17,20,21,22)],
            est = data4$beta_IVW,
            lower = data4$lowerCI95, 
            upper = data4$upperCI95,
            sizes = 0.5,
            ci_column = 2,
            ref_line = 0,
            #title = myXlab,
            #xlab = myXlab,
            xlim = c(-1,7),
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/SupFigs/ForestPlot_2v1_sampleMR.png")
png(filename = filename,width = 1900, height = 700, res=200)
plot(p2)
dev.off()

#' ## All SNPs vs nominal significant SNPs ####
#' 
MVMR2 = copy(MVMR)
MVMR2 = MVMR2[outcome == "UKB_BW",]
MVMR2 = MVMR2[flag == "0",]
MVMR2 = MVMR2[exposure == "logefwcomb",]
MVMR2 = MVMR2[setting == "multivariate",]

setorder(MVMR2,flag)
MVMR2[,rank := 2]
dummy = data.table(exposure_type = unique(MVMR2$exposure_type),rank=1)
data4 = rbind(MVMR2,dummy, fill=T)
data4[,subgroup := paste0("   ", flag)]
data4[is.na(threshold),subgroup := exposure_type]

data4[,lowerCI95 := beta_IVW-1.96*SE_IVW]
data4[,upperCI95 := beta_IVW+1.96*SE_IVW]
data4$` ` <- paste(rep(" ", 50), collapse = " ")
data4$`Estimate \n[95% CI]` <- ifelse(is.na(data4$SE_IVW), "",
                                      sprintf("%.2f [%.2f, %.2f]",
                                              data4$beta_IVW, data4$lowerCI95, data4$upperCI95))

setorder(data4,exposure_type,rank)
data4[exposure_type=="var" & is.na(setting),subgroup:= "variability"]
data4[grepl("all",threshold), subgroup := gsub("0","main (all SNPs)",subgroup)]
data4[grepl("nominal",threshold), subgroup := gsub("0","sensitivity (nom. sig. SNPs)",subgroup)]

dummy = c("white","#F2AA84", "#FBE3D6",
          "white","#47D45A", "#C2F1C8", 
          "white","#61CBF4", "#CAEEFB" )
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab =  "increase in BW (95% CI) per 1-SD increment EFW"

data4[,condF := round(condFstat,1)]
data4[,condF := as.character(condF)]
data4[is.na(setting),condF := ""]
setnames(data4,"condF","cond. \nF-stat")
setnames(data4,"subgroup", "Exposure type \n   SNP selection")

p2<- forest(data4[,c(17,20,21,22)],
            est = data4$beta_IVW,
            lower = data4$lowerCI95, 
            upper = data4$upperCI95,
            sizes = 0.5,
            ci_column = 2,
            ref_line = 0,
            #title = myXlab,
            #xlab = myXlab,
            xlim = c(-1,6),
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/SupFigs/ForestPlot_allVnom_SNPs.png")
png(filename = filename,width = 1900, height = 700, res=200)
plot(p2)
dev.off()

#' ## log vs Z-score ####
MVMR2 = copy(MVMR)
MVMR2 = MVMR2[outcome == "UKB_BW",]
MVMR2 = MVMR2[threshold == "all_SNPs",]
MVMR2 = MVMR2[flag %in%  c("0","1A")]
MVMR2 = MVMR2[exposure %in% c("logefwcomb","efwcombZv2"),]
MVMR2[,rank := 2]

dummy = data.table(exposure_type = unique(MVMR2$exposure_type),rank=1)
data4 = rbind(MVMR2,dummy, fill=T)
data4[,subgroup := paste0("   ", flag, " - ",setting)]
data4[,subgroup := gsub("0", "main", subgroup)]
data4[,subgroup := gsub("1A", "no slope", subgroup)]
data4[,subgroup := gsub("1B", "no variabiltity", subgroup)]
data4[,subgroup := gsub("2", "no statins", subgroup)]
data4[,subgroup := gsub("multivariate", "MVMR", subgroup)]
data4[,subgroup := gsub("univariate", "     MR", subgroup)]
data4[is.na(threshold),subgroup := exposure_type]

data4[,lowerCI95 := beta_IVW-1.96*SE_IVW]
data4[,upperCI95 := beta_IVW+1.96*SE_IVW]
data4$ES <- paste(rep(" ", 15), collapse = " ")
data4$CI <- ifelse(is.na(data4$SE_IVW), "",sprintf("%.2f [%.2f, %.2f]",data4$beta_IVW, data4$lowerCI95, data4$upperCI95))

setorder(data4,exposure_type,rank)
data4[exposure_type=="var" & is.na(setting),subgroup:= "variability"]
data4[,condF := round(condFstat,1)]
data4[,condF := as.character(condF)]
data4[is.na(setting),condF := ""]

data5 = copy(data4)
data5 = data5[is.na(exposure) | exposure == "logefwcomb",c(17,20,21,22,7,18,19)]
data6 = copy(data4)
data6 = data6[is.na(exposure) | exposure == "efwcombZv2",c(17,20,21,22,7,18,19)]

data7 = cbind(data5,data6[,2:7])
names(data7)
names(data7) = c("Exposure type \n   Setting", 
                 "EFW \nlog-transformed", "Estimate \n[95% CI]","(cond.) \nF-stat","beta_log","lCI_log","uCI_log",
                 "EFW \nZ-score", "Estimate \n[95% CI]","(cond.) \nF-stat","beta_Z","lCI_Z","uCI_Z")

dummy = c("white","#F2AA84", rep("#FBE3D6",3),
          "white","#47D45A", "#C2F1C8",
          "white","#61CBF4", rep("#CAEEFB",3))
tm1<- forest_theme(base_size = 10,
                   core=list(bg_params=list(fill = dummy)))

myXlab =  "increase in BW (95% CI) per 1-SD increment in EFW levels, weekly increase, and variability"

p2<- forest(data7[,c(1,2,3,4,8,9,10)],
            est = list(data7$beta_log,
                       data7$beta_Z),
            lower = list(data7$lCI_log,
                         data7$lCI_Z), 
            upper = list(data7$uCI_log,
                         data7$uCI_Z),
            ci_column = c(2,5),
            sizes = 0.5,
            ref_line = 0,
            ticks_at = list(c(-1, 1, 3,5), c(0,0.2,0.4)),
            #title = myXlab,
            #xlab = c("",myXlab),
            #xlim = c(-0.5,6),
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/SupFigs/ForestPlot_differentExposures.png")
png(filename = filename,width = 1600, height = 750, res=200)
plot(p2)
dev.off()


#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))
