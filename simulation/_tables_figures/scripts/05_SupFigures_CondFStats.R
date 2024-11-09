#' ---
#' title: "Supplemental simulation figure: conditional F-statistics"
#' subtitle: "Supplemental Figures"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output: html_document
#' ---
#'
#' # Introduction ####
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../../SourceFile_HPC.R")

outdir_results = "../results/_figures/SupFigures/"
if(dir.exists(outdir_results)==F){
  dir.create(outdir_results)
  message("Created results folder ",outdir_results, " for main figures ")
}else{
  message("Using pre-existing results folder ",outdir_results, " for main figures ")
}

#' # Get data ####
#' ***
load("../results/SupTabs.RData")

myTab = copy(tab2)
myTab = myTab[Sim_Y %in% c("CM2")]
myTab = myTab[outcome == "Y8"]

#' # Get plots ####
#' ***
dumTab1 = copy(myTab)
dumTab1[,type := "mean"]

dumTab2 = copy(dumTab1)
dumTab2[,condFStats_median_mean := condFStats_median_slope]
dumTab2[,condFStats_1stQ_mean := condFStats_1stQ_slope]
dumTab2[,condFStats_3rdQ_mean := condFStats_3rdQ_slope]
dumTab2[,type := "slope"]

dumTab3 = copy(dumTab1)
dumTab3[,condFStats_median_mean := condFStats_median_var]
dumTab3[,condFStats_1stQ_mean := condFStats_1stQ_var]
dumTab3[,condFStats_3rdQ_mean := condFStats_3rdQ_var]
dumTab3[,type := "variability"]

dumTab4 = rbind(dumTab1,dumTab2,dumTab3)  
dumTab4 = dumTab4[!is.na(condFStats_median_mean),]
dumTab4[,dumID := paste(Sim_X,Sim_Y, sep=" - ")]

dumTab4[,check := gsub("main_","",check)]
dumTab4[,check := gsub("sens_","",check)]
dumTab4[,check := gsub("_"," - ",check)]

plot_X12 = ggplot(dumTab4[Sim_X == "X12"], aes(x=check, y=condFStats_median_mean, col = type,shape = type)) +
  geom_hline(yintercept = 0,color="grey") +
  geom_hline(yintercept = 10,color="black",linetype = "dashed") +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=condFStats_1stQ_mean, 
                    ymax=condFStats_3rdQ_mean), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme(legend.position = "none") +
  xlab("Scenarios") + ylab("Conditional F-Statistics") 
plot_X12

plot_X123 = ggplot(dumTab4[Sim_X == "X123"], aes(x=check, y=condFStats_median_mean, col = type,shape = type)) +
  geom_hline(yintercept = 0,color="grey") +
  geom_hline(yintercept = 10,color="black",linetype = "dashed") +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=condFStats_1stQ_mean, 
                    ymax=condFStats_3rdQ_mean), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme(legend.position = "none") +
  xlab("Scenarios") + ylab("Conditional F-Statistics") 
plot_X123

plot_X13 = ggplot(dumTab4[Sim_X == "X13"], aes(x=check, y=condFStats_median_mean, col = type,shape = type)) +
  geom_hline(yintercept = 0,color="grey") +
  geom_hline(yintercept = 10,color="black",linetype = "dashed") +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=condFStats_1stQ_mean, 
                    ymax=condFStats_3rdQ_mean), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme(legend.position = "none") +
  xlab("Scenarios") + ylab("Conditional F-Statistics") 
plot_X13

dumTab5 = copy(dumTab4)
dumTab5 = dumTab5[check %in% c("A","B","A - noSlope","B - noSlope"),]

dumTab5[,type2:= "main"]
dumTab5[grepl("noSlope",check),type2:= "no slope"]

dumTab5[,flag2:= "A"]
dumTab5[grepl("B",check),flag2:= "B"]

plot_main = ggplot(dumTab5, aes(x=flag2, y=condFStats_median_mean, col = type2)) +
  facet_grid2(type~ Sim_X,scales = "free",independent = "y") +
  geom_hline(yintercept = 0,color="grey") +
  geom_hline(yintercept = 10,color="black",linetype = "dashed") +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=condFStats_1stQ_mean, 
                    ymax=condFStats_3rdQ_mean), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  #theme(legend.position = "none") +
  xlab("scenarios") + ylab("conditional F-Statistics") + labs(col = "GAMLSS \nmodel")
plot_main

filename = paste0(outdir_results,"/CondFStats_X123.png")
png(filename = filename,width = 3200, height = 1600, res=200)
print(plot_X123)
dev.off()

filename = paste0(outdir_results,"/CondFStats_X13.png")
png(filename = filename,width = 3200, height = 1600, res=200)
print(plot_X13)
dev.off()

filename = paste0(outdir_results,"/CondFStats_X12.png")
png(filename = filename,width = 3200, height = 1600, res=200)
print(plot_X12)
dev.off()

filename = paste0(outdir_results,"/CondFStats_main.png")
png(filename = filename,width = 2400, height = 1400, res=200)
print(plot_main)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
