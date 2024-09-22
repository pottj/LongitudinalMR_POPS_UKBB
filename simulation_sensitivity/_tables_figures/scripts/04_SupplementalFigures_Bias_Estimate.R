#' ---
#' title: "Supplemental Figures for simulation - bias"
#' subtitle: "Longitudinal MVMR"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output: html_document
#' ---
#'
#' # Introduction ####
#' ***
#' I want panel plots for the main simulation picturing bias and corrected estiamtes next to each other. Let's see if this works. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../../SourceFile_HPC.R")

#' # Load data #### 
#' ***
load("../results/_tables/Simulation_complete_main_v1.RData")

#' # Get Bias table ####
#' ***
dumTab = copy(myTab)
dumTab[,type := "mean"]

dumTab2 = copy(dumTab)
dumTab2[,bias_X1 := bias_X2]
dumTab2[,bias_SE_X1 := bias_SE_X2]
dumTab2[,type := "slope"]

dumTab3 = copy(dumTab)
dumTab3[,bias_X1 := bias_X3]
dumTab3[,bias_SE_X1 := bias_SE_X3]
dumTab3[,type := "variability"]

dumTab4 = rbind(dumTab,dumTab2,dumTab3)
dumTab4 = dumTab4[!is.na(bias_X1),]

# reduce to relevant columns
myNames = c("Sim_X", "Sim_Y", "outcome", "type", "bias_X1", "bias_SE_X1")
colsOut = setdiff(names(dumTab4),myNames)
dumTab4[,get("colsOut"):=NULL]
setcolorder(dumTab4,myNames)
setnames(dumTab4,"bias_X1","value")
setnames(dumTab4,"bias_SE_X1","SE")
dumTab4[,type := paste0("bias - ",type)]

biasData = copy(dumTab4)

#' # Get corrected estimate table ####
#' ***
dumTab = copy(myTab)
dumTab[,type := "mean"]

dumTab2 = copy(dumTab)
dumTab2[,mean_betaIVW_X1 := mean_betaIVW_X2]
dumTab2[,empSE_X1 := empSE_X2]
dumTab2[,type := "slope"]

dumTab3 = copy(dumTab)
dumTab3[,mean_betaIVW_X1 := mean_betaIVW_X3]
dumTab3[,empSE_X1 := empSE_X3]
dumTab3[,type := "variability"]

dumTab4 = rbind(dumTab,dumTab2,dumTab3)
dumTab4 = dumTab4[!is.na(mean_betaIVW_X1),]

# reduce to relevant columns
myNames = c("Sim_X", "Sim_Y", "outcome", "type", "mean_betaIVW_X1", "empSE_X1")
colsOut = setdiff(names(dumTab4),myNames)
dumTab4[,get("colsOut"):=NULL]
setcolorder(dumTab4,myNames)
setnames(dumTab4,"mean_betaIVW_X1","value")
setnames(dumTab4,"empSE_X1","SE")
dumTab4[,type := paste0("estimate - ",type)]

estimateData = copy(dumTab4)

#' # Get Plots ####
#' ***
plotData = rbind(biasData,estimateData)

data_hlines = data.frame(type = c(rep("estimate - mean",4),"estimate - slope","estimate - variability"),
                         mylines = c(0.3,1.2,-1.2,-0.3,0.3,1))

plot12 = ggplot(plotData[Sim_X == "X12",], aes(x=Sim_Y, y=value, color = outcome)) +
  facet_wrap(~ type,scales = "free_y") +
  geom_hline(yintercept = 0,color="grey") +
  geom_hline(data = data_hlines, col="darkred", linetype="dotted", 
             linewidth=1.3, aes(yintercept = mylines)) +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=value-1.96*SE, ymax=value+1.96*SE), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  xlab("Causal model") + ylab("") +
  labs(color = "Outcome")
plot12

filename = paste0("../results/SupFig_Bias_Estimate_X12.png")
png(filename = filename,width = 2800, height = 1600, res=200)
print(plot12)
dev.off()

data_hlines2 = data.frame(type = c(rep("estimate - mean",5)),
                          mylines = c(-0.9,-0.6,0.3,-2.1,-1.2))

plot13 = ggplot(plotData[Sim_X == "X13",], aes(x=Sim_Y, y=value, color = outcome)) +
  facet_wrap(~ type,scales = "free_y") +
  geom_hline(yintercept = 0,color="grey") +
  geom_hline(data = data_hlines, col="darkred", linetype="dotted", 
             linewidth=1.3, aes(yintercept = mylines)) +
  geom_hline(data = data_hlines2, col="black", linetype="dashed", 
             linewidth=0.5, aes(yintercept = mylines)) +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=value-1.96*SE, ymax=value+1.96*SE), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  xlab("Causal model") + ylab("") +
  labs(color = "Outcome")
plot13

filename = paste0("../results/SupFig_Bias_Estimate_X13.png")
png(filename = filename,width = 2800, height = 1600, res=200)
print(plot13)
dev.off()

plot123 = ggplot(plotData[Sim_X == "X123",], aes(x=Sim_Y, y=value, color = outcome)) +
  facet_wrap(~ type,scales = "free_y") +
  geom_hline(yintercept = 0,color="grey") +
  geom_hline(data = data_hlines, col="darkred", linetype="dotted", 
             linewidth=1.3, aes(yintercept = mylines)) +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=value-1.96*SE, ymax=value+1.96*SE), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  xlab("Causal model") + ylab("") +
  labs(color = "Outcome")
plot123

filename = paste0("../results/SupFig_Bias_Estimate_X123.png")
png(filename = filename,width = 2800, height = 1600, res=200)
print(plot123)
dev.off()

#' # Session Info ####
#' ***

sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")


