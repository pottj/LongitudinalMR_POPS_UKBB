#' ---
#' title: "Get Supplemental Figures for real data (forest plots)"
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
#' **Supplemental Figures: Forest Plots**
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
myFiles = list.files(path = "../results/", pattern = "03_MVMR")
mySettings = gsub(".RData","",myFiles)
mySettings = gsub("03_MVMR_0._","",mySettings)
mySettings
mySettings2 = c("0","1A","1B","2A","2B","2C","3A","3B","3C")

dumTab2 = foreach(i = 1:length(myFiles))%do%{
  #i=1
  load(paste0("../results/",myFiles[i]))
  MVMR = copy(MVMR_results)
  #MVMR = MVMR[outcome == "Aragam"]
  MVMR = MVMR[setting == "multivariate"]
  MVMR = MVMR[threshold == "all_SNPs"]
  MVMR[,flag := mySettings2[i]]
  MVMR
}

MVMR = rbindlist(dumTab2)

#' # Plotting ####
#' ***
#' ## Different Outcome data
#' 
MVMR2 = copy(MVMR)
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
data4[grepl("Aragam",outcome), subgroup := gsub("0","main (Aragam outcome data)",subgroup)]
data4[grepl("UKB",outcome), subgroup := gsub("0","sensitivity (UKB outcome data)",subgroup)]
setorder(data4,exposure_type,rank,subgroup)

dummy = c("white","#F2AA84", "#FBE3D6",
          "white","#47D45A", "#C2F1C8", 
          "white","#61CBF4", "#CAEEFB" )
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab =  "logOR for CAD risk (95% CI) per 1-SD increment in TC levels, yearly increase, and variability"

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
            xlim = c(-1,1),
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/SupFigs/ForestPlot_2v1_sampleMR.png")
png(filename = filename,width = 1900, height = 700, res=200)
plot(p2)
dev.off()

#' ## Different exposure data
MVMR2 = copy(MVMR)
MVMR2 = MVMR2[flag %in% c("0","2A","2B","2C")]
MVMR2 = MVMR2[outcome == "Aragam"]
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
data4[grepl("0",subgroup), subgroup := gsub("0","main analysis",subgroup)]
data4[grepl("2A",subgroup), subgroup := gsub("2A","2A - no statins",subgroup)]
data4[grepl("2B",subgroup), subgroup := gsub("2B","2B - after BL",subgroup)]
data4[grepl("2C",subgroup), subgroup := gsub("2C","2C - before BL",subgroup)]

dummy = c("white","#F2AA84", "#FBE3D6", "#FBE3D6", "#FBE3D6",
          "white","#47D45A", "#C2F1C8", "#C2F1C8", "#C2F1C8", 
          "white","#61CBF4", "#CAEEFB", "#CAEEFB", "#CAEEFB")
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab =  "logOR for CAD risk (95% CI) per 1-SD increment in TC levels, yearly increase, and variability"

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
            xlim = c(-1,1),
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/SupFigs/ForestPlot_SampleSetting_Aragam.png")
png(filename = filename,width = 1900, height = 1000, res=200)
plot(p2)
dev.off()

#' ## Different exposure data - 1 sample
MVMR2 = copy(MVMR)
MVMR2 = MVMR2[flag %in% c("0","2A","2B","2C")]
MVMR2 = MVMR2[outcome != "Aragam"]
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
data4[grepl("0",subgroup), subgroup := gsub("0","main analysis",subgroup)]
data4[grepl("2A",subgroup), subgroup := gsub("2A","2A - no statins",subgroup)]
data4[grepl("2B",subgroup), subgroup := gsub("2B","2B - after BL",subgroup)]
data4[grepl("2C",subgroup), subgroup := gsub("2C","2C - before BL",subgroup)]

dummy = c("white","#F2AA84", "#FBE3D6", "#FBE3D6", "#FBE3D6",
          "white","#47D45A", "#C2F1C8", "#C2F1C8", "#C2F1C8", 
          "white","#61CBF4", "#CAEEFB", "#CAEEFB", "#CAEEFB")
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab =  "logOR for CAD risk (95% CI) per 1-SD increment in TC levels, yearly increase, and variability"

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
            xlim = c(-1,1.5),
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/SupFigs/ForestPlot_SampleSetting_UKB.png")
png(filename = filename,width = 1900, height = 1000, res=200)
plot(p2)
dev.off()

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))

