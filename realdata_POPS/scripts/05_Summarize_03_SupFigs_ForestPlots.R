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
#' - main: 2-sample (UKB) vs 1-sample (POPS)
#' - main: all SNPs vs SNP selection
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
  #MVMR = MVMR[outcome == "UKB_BW"]
  MVMR = MVMR[setting == "multivariate"]
  #MVMR = MVMR[threshold == "all_SNPs"]
  MVMR[,flag := mySettings2[i]]
  MVMR
}

MVMR = rbindlist(dumTab2)

#' # Plotting ####
#' ***
#' ## 2-sample vs 1-sample ####
#' 
MVMR2 = copy(MVMR)
MVMR2 = MVMR2[outcome != "EGG_BW",]
MVMR2 = MVMR2[threshold == "all_SNPs",]
MVMR2 = MVMR2[flag %in% c("0")]

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
data4[grepl("UKB",outcome), subgroup := gsub("0","main (UKB outcome data)",subgroup)]
data4[grepl("POPS",outcome), subgroup := gsub("0","sensitivity (POPS outcome data)",subgroup)]
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
setnames(data4,"subgroup", "exposure type \n   MR approach")

p2<- forest(data4[,c(17,20,21,22)],
            est = data4$beta_IVW,
            lower = data4$lowerCI95, 
            upper = data4$upperCI95,
            sizes = 0.5,
            ci_column = 2,
            ref_line = 0,
            #title = myXlab,
            xlab = myXlab,
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
setnames(data4,"subgroup", "exposure type \n   MR approach")

p2<- forest(data4[,c(17,20,21,22)],
            est = data4$beta_IVW,
            lower = data4$lowerCI95, 
            upper = data4$upperCI95,
            sizes = 0.5,
            ci_column = 2,
            ref_line = 0,
            #title = myXlab,
            xlab = myXlab,
            xlim = c(-1,6),
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/SupFigs/ForestPlot_allVnom_SNPs.png")
png(filename = filename,width = 1900, height = 700, res=200)
plot(p2)
dev.off()

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))

