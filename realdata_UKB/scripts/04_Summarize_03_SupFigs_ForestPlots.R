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
mySettings2 = c("0","1B","1A","2","3")

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
data4[grepl("Aragam",outcome), subgroup := gsub("0","main (Aragam)",subgroup)]
data4[grepl("UKB",outcome), subgroup := gsub("0","sensitivity (UKB)",subgroup)]
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
            #xlim = c(-1,1),
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/SupFigs/ForestPlot_2v1_sampleMR.png")
png(filename = filename,width = 1700, height = 650, res=200)
plot(p2)
dev.off()

#' ## SNP selection ####
#' 
MVMR2 = copy(MVMR)
MVMR2 = MVMR2[outcome == "Aragam",]
MVMR2 = MVMR2[setting == "multivariate",]
MVMR2 = MVMR2[(flag=="3" & threshold=="all_SNPs") | flag=="0"]
MVMR2[flag=="0" & threshold == "gw_SNPs", ID := "sens_SNPset_gw"]
MVMR2[,rank := 2]

dummy = data.table(exposure_type = unique(MVMR2$exposure_type),rank=1)
data4 = rbind(MVMR2,dummy, fill=T)
data4[,subgroup := paste0("   ", ID)]
data4[is.na(threshold),subgroup := exposure_type]

data4[,lowerCI95 := beta_IVW-1.96*SE_IVW]
data4[,upperCI95 := beta_IVW+1.96*SE_IVW]
data4$` ` <- paste(rep(" ", 50), collapse = " ")
data4$`Estimate \n[95% CI]` <- ifelse(is.na(data4$SE_IVW), "",
                                      sprintf("%.2f [%.2f, %.2f]",
                                              data4$beta_IVW, data4$lowerCI95, data4$upperCI95))

setorder(data4,exposure_type,rank)
data4[exposure_type=="var" & is.na(setting),subgroup:= "variability"]
data4[, subgroup := gsub("main","main analysis",subgroup)]
data4[, subgroup := gsub("sens_SNPset_gw","  only gw. sig. SNPs",subgroup)]
data4[, subgroup := gsub("sens_SNPset_distinct","  distinct sets for M and V",subgroup)]
data4[, subgroup := gsub("sens_SNPset_top20","  top20 for M, S, and V",subgroup)]
data4[, subgroup := gsub("sens_SNPset_pleiotropy","  only pleiotropic SNPs",subgroup)]
data4[, subgroup := gsub("sens_SNPset_bioknown","  biological plausible SNPs",subgroup)]
data4[, subgroup := gsub("sens_SNPset_sexIA","  SNP with sig. sex IA",subgroup)]

dummy = c("white","#F2AA84", rep("#FBE3D6", 6),
          "white","#47D45A", rep("#C2F1C8", 6), 
          "white","#61CBF4", rep("#CAEEFB", 6))
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab =  "logOR for CAD risk (95% CI) per 1-SD increment in TC levels, yearly increase, and variability"

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
            xlim = c(-1.5,2),
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/SupFigs/ForestPlot_SNPsetting_Aragam.png")
png(filename = filename,width = 1900, height = 1400, res=200)
plot(p2)
dev.off()

#' ## MVMR vs MR ####
MVMR2 = copy(MVMR)
MVMR2 = MVMR2[outcome == "Aragam",]
MVMR2 = MVMR2[threshold == "all_SNPs",]
MVMR2 = MVMR2[flag %in%  c("0","1A","2")]
MVMR2[,rank := 2]

dummy = data.table(exposure_type = unique(MVMR2$exposure_type),rank=1)
data4 = rbind(MVMR2,dummy, fill=T)
data4[,subgroup := paste0("   ", flag, " - ",setting)]
data4[,subgroup := gsub("0", "main", subgroup)]
data4[,subgroup := gsub("1A", "no slope", subgroup)]
data4[,subgroup := gsub("2", "no statins", subgroup)]
data4[,subgroup := gsub("multivariate", "MVMR", subgroup)]
data4[,subgroup := gsub("univariate", "     MR", subgroup)]
data4[is.na(threshold),subgroup := exposure_type]

data4[,lowerCI95 := beta_IVW-1.96*SE_IVW]
data4[,upperCI95 := beta_IVW+1.96*SE_IVW]
data4$ES <- paste(rep(" ", 25), collapse = " ")
data4$CI <- ifelse(is.na(data4$SE_IVW), "",sprintf("%.2f [%.2f, %.2f]",data4$beta_IVW, data4$lowerCI95, data4$upperCI95))

setorder(data4,exposure_type,rank)
data4[exposure_type=="var" & is.na(setting),subgroup:= "variability"]
data4[,condF := round(condFstat,1)]
data4[,condF := as.character(condF)]
data4[is.na(setting),condF := ""]

data5 = copy(data4)
data5 = data5[is.na(exposure) | setting == "multivariate",c(17,20,21,22,7,18,19)]
data6 = copy(data4)
data6 = data6[is.na(exposure) | setting == "univariate",c(17,20,21,22,7,18,19)]

data7 = cbind(data5,data6[,2:7])
data7[,subgroup := gsub(" - MVMR","",subgroup)]
names(data7)
names(data7) = c("Exposure type \n   Setting", 
                 "\nMVMR", "MVMR Estimate \n[95% CI]","cond. \nF-stat","beta_MVMR","lCI_MVMR","uCI_MVMR",
                 "\nMR", "MR Estimate \n[95% CI]","\nF-stat","beta_MR","lCI_MR","uCI_MR")

dummy = c("white","#F2AA84", rep("#FBE3D6",2),
          "white","#47D45A", "#C2F1C8",
          "white","#61CBF4", rep("#CAEEFB",2))
tm1<- forest_theme(base_size = 10,
                   core=list(bg_params=list(fill = dummy)))

p2<- forest(data7[,c(1,2,3,4,8,9,10)],
            est = list(data7$beta_MVMR,
                       data7$beta_MR),
            lower = list(data7$lCI_MVMR,
                         data7$lCI_MR), 
            upper = list(data7$uCI_MVMR,
                         data7$uCI_MR),
            ci_column = c(2,5),
            sizes = 0.5,
            ref_line = 0,
            ticks_at = list(c(-0.5, -0.25, 0, 0.25,0.5), c(-0.25,0,0.5,1)),
            #title = myXlab,
            #xlab = c("",myXlab),
            #xlim = c(-0.5,6),
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/SupFigs/ForestPlot_MVMRvMR.png")
png(filename = filename,width = 1900, height = 700, res=200)
plot(p2)
dev.off()


#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))

