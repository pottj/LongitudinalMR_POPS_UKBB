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

source("../../SourceFile.R")

#' # Load data ####
#' ***
load(paste0(UKB_phenotypes_filtered,"/01_Prep_01_BL_FU_GP_merged_filtered.RData"))

#' # Get plot data ####
#' ***
plotData = copy(myTab_long)
names(plotData)

plotData[,table(sex)]
plotData[sex==2,sex2 := "female"]
plotData[sex==1,sex2 := "male"]

#' # Get Plot 1: all samples ####
#' ***
ggp1  = ggplot(plotData, aes(x=age, y=TC, group=ID,
                             col=sex2,fill=sex2,
                             shape = as.factor(statin))) +
  geom_hline(yintercept = 8,linetype="dashed", linewidth=0.5)+
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

filename = paste0("../results/_figures/SupFigs/Trajectory_TC_MAIN.png")
png(filename = filename,width=1500,height=1200,res = 200)
plot(ggp1)
dev.off()

#' # Get Plot 2: 10 random samples ####
#' ***
set.seed(2023)
myIDs = sample(unique(plotData$ID),size=15)

ggp2 = ggplot(plotData[ID %in% myIDs], 
                 aes(x=age, y=TC, group=ID,shape=as.factor(statin))) +
  geom_hline(yintercept = 8,linetype="dotted", linewidth=0.5)+
  geom_line(aes(col=as.factor(ID))) + 
  geom_point(aes(colour = as.factor(ID)))+
  labs(x="Age (in years)",
       y="Total cholesterol (in mmol/l)", 
       color="IDs",
       shape="Statin\ntreatment",
       title = paste0("Trajectories of TC levels (example with 15 random samples)")) +
  scale_shape_manual(values=c(21,24),
                     labels = c("no","yes"))+
  # scale_colour_manual(values = c("darkred","steelblue"),
  #                     labels = c("women","men"))+
  theme_classic() 

filename = paste0("../results/_figures/SupFigs/Trajectory_TC_MAIN_randomSamples.png")
png(filename = filename,width = 1800, height = 1000, res=200)
plot(ggp2)
dev.off()

#' # Get Plot 3: sensitivity subset ####
#' ***
ggp3  = ggplot(plotData[sens1==T,], aes(x=age, y=TC, group=ID,
                             col=sex2,fill=sex2)) +
  geom_hline(yintercept = 8,linetype="dotted", linewidth=0.5)+
  geom_line(aes(alpha=0.01)) + 
  geom_point(shape = 21,color="black")+
  labs(x="Age (years)", y="TC (mmol/L)") +
  scale_fill_manual(values = c("darkred","steelblue"),
                    labels = c("female", "male"))+
  scale_colour_manual(values = c("darkred","steelblue"),
                      labels = c("female", "male"))+
  theme_classic() +
  theme(legend.position = "none")

filename = paste0("../results/_figures/SupFigs/Trajectory_TC_SENS_sampleSet.png")
png(filename = filename,width=1500,height=1200,res = 200)
plot(ggp3)
dev.off()

#' # Get Plot 2: 10 random samples ####
#' ***
set.seed(2023)
myIDs = sample(unique(plotData[sens1==T,ID]),size=15)

ggp2 = ggplot(plotData[ID %in% myIDs & sens1==T,], 
              aes(x=age, y=TC, group=ID,shape=as.factor(statin))) +
  geom_hline(yintercept = 8,linetype="dotted", linewidth=0.5)+
  geom_line(aes(col=as.factor(ID))) + 
  geom_point(aes(colour = as.factor(ID)))+
  labs(x="Age (in years)",
       y="Total cholesterol (in mmol/l)", 
       color="IDs",
       shape="Statin\ntreatment",
       title = paste0("Trajectories of TC levels (example with 15 random samples)")) +
  scale_shape_manual(values=c(21,24),
                     labels = c("no","yes"))+
  # scale_colour_manual(values = c("darkred","steelblue"),
  #                     labels = c("women","men"))+
  theme_classic() 

filename = paste0("../results/_figures/SupFigs/Trajectory_TC_SENS_randomSamples.png")
png(filename = filename,width = 1800, height = 1000, res=200)
plot(ggp2)
dev.off()

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))


