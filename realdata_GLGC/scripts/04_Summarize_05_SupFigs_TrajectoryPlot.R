#' ---
#' title: "Get Supplemental Figures for real data (trajectory plot)"
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
#' **Supplemental Figures: Trajectory Plot**
#' 
#' - MAIN (all samples): colored by sex
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../../SourceFile_HPC.R")

#' # Load data ####
#' ***
load(paste0(UKB_phenotypes_filtered,"/01_Prep_01_UKB_GP_TC_GLGC.RData"))

#' # Get plot data ####
#' ***
plotData = copy(myTab6)
names(plotData)
matched = match(plotData$BSU_ID,myTab7$ID)
plotData[,sex := myTab7[matched,sex]]
plotData[sex==0,sex2 := "female"]
plotData[sex==1,sex2 := "male"]

#' # Get Plot ####
#' ***
ggp2  = ggplot(plotData, aes(x=exposure_age, y=exposure_value, group=BSU_ID,
                             col=sex2,fill=sex2,
                             shape = as.factor(lipLowMed))) +
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  geom_line(aes(alpha=0.01)) + 
  geom_point(color="black")+
  labs(x="Age (years)", y="TC (mmol/L)") +
  scale_shape_manual(values = c(21,24),
                     labels = c("no treatment","statin treatment"))+
  scale_fill_manual(values = c("darkred","steelblue"),
                      labels = c("female", "male"))+
  scale_colour_manual(values = c("darkred","steelblue"),
                      labels = c("female", "male"))+
  theme_classic() +
  theme(legend.position = "none")
ggp2

filename = paste0("../results/_figures/SupFigs/Trajectory_TC_MAIN.png")
png(filename = filename,width=1500,height=1200,res = 200)
plot(ggp2)
dev.off()

#' # Load data ####
#' ***
load(paste0(UKB_phenotypes_filtered,"/01_Prep_05_UKB_GP_TC_GLGC_sens.RData"))

#' # Get plot data ####
#' ***
plotData = copy(myTab6)
names(plotData)
matched = match(plotData$BSU_ID,myTab7$ID)
plotData[,sex := myTab7[matched,sex]]
plotData[sex==0,sex2 := "female"]
plotData[sex==1,sex2 := "male"]

#' # Get Plot ####
#' ***
ggp2  = ggplot(plotData, aes(x=exposure_age, y=exposure_value, group=BSU_ID,
                             col=sex2,fill=sex2)) +
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  geom_line(aes(alpha=0.01)) + 
  geom_point(shape = 21,color="black")+
  labs(x="Age (years)", y="TC (mmol/L)") +
  scale_fill_manual(values = c("darkred","steelblue"),
                    labels = c("female", "male"))+
  scale_colour_manual(values = c("darkred","steelblue"),
                      labels = c("female", "male"))+
  theme_classic() +
  theme(legend.position = "none")
ggp2

filename = paste0("../results/_figures/SupFigs/Trajectory_TC_SENS_sampleSet.png")
png(filename = filename,width=1500,height=1200,res = 200)
plot(ggp2)
dev.off()

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))


