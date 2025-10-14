#' ---
#' title: "Get Forest Plots"
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
mySettings2 = c("0","2","1A","1B")

dumTab2 = foreach(i = 1:length(myFiles))%do%{
  #i=1
  load(paste0("../results/",myFiles[i]))
  MVMR = copy(MVMR_results)
  MVMR = MVMR[outcome == "UKB_BW"]
  MVMR = MVMR[setting == "multivariate"]
  MVMR = MVMR[threshold == "all_SNPs"]
  MVMR[,flag := mySettings2[i]]
  MVMR
}

MVMR = rbindlist(dumTab2)

#' # Plotting ####
#' ***
#' ## Main and no slope/var tests ####
#' 
MVMR2 = copy(MVMR)
MVMR2 = MVMR2[flag %in% c("0","1A","1B")]
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
data4[grepl("1A",subgroup), subgroup := gsub("1A","1A - no variability",subgroup)]
data4[grepl("1B",subgroup), subgroup := gsub("1B","1B - no slope",subgroup)]

dummy = c("white","#F2AA84", "#FBE3D6", "#FBE3D6",
          "white","#47D45A", "#C2F1C8", 
          "white","#61CBF4","#CAEEFB")
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

filename = paste0("../results/_figures/MainFigures/EFW_BW_main_vs_GAMLSSmod.png")
png(filename = filename,width = 1700, height = 750, res=200)
plot(p2)
dev.off()

#' ## Main and different sample settings ####
MVMR2 = copy(MVMR)
MVMR2 = MVMR2[flag %in% c("0","2")]
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
data4[grepl("2",subgroup), subgroup := gsub("2","2 - GBR subset",subgroup)]

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

filename = paste0("../results/_figures/MainFigures/EFW_BW_main_vs_SampleSetting.png")
png(filename = filename,width = 1900, height = 1000, res=200)
plot(p2)
dev.off()

#' ## All in one ####
MVMR2 = copy(MVMR)
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
data4[grepl("0",subgroup), subgroup := gsub("0","main analysis",subgroup)]
data4[grepl("1A",subgroup), subgroup := gsub("1A","1A - no variability",subgroup)]
data4[grepl("1B",subgroup), subgroup := gsub("1B","1B - no slope",subgroup)]
data4[grepl("2",subgroup), subgroup := gsub("2","2 - GBR subset",subgroup)]

dummy = c("white","#F2AA84", "#FBE3D6", "#FBE3D6", "#FBE3D6", 
          "white","#47D45A", "#C2F1C8", "#C2F1C8", 
          "white","#61CBF4", "#CAEEFB", "#CAEEFB" )
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

filename = paste0("../results/_figures/MainFigures/EFW_BW_main_vs_additionalSettings.png")
png(filename = filename,width = 1900, height = 1000, res=200)
plot(p2)
dev.off()


#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))
