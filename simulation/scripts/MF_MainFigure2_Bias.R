#' ---
#' title: "Evaluation POPS simulation - main figure 2"
#' subtitle: "Summary of all simulations"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output: html_document
#' ---
#'
#' # Introduction ####
#' ***
#' For the paper, I want one figure for the bias
#' 
#' - only use quadratic growth
#' - only use mean and slope
#' - only use Y2, Y3, and Y4 & their binary counterparts
#' 
#' # Initialize ####
#' ***

rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_HPC.R")

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("202.-","24-",tag)
tag = gsub("-","",tag)

#' # Load data ####
#' ***
load("../results/_tables/Simulation_complete.RData")

#' # Get figure ####
#' ***
#' Filter data 
dumTab = copy(myTab)
dumTab = dumTab[Sim_growth=="quad"]
dumTab = dumTab[!grepl("Y1",outcome),]

dumTab[,type := "mean"]

dumTab2 = copy(dumTab)
dumTab2[,bias_X1 := bias_X2]
dumTab2[,bias_SE_X1 := bias_SE_X2]
dumTab2[,type := "slope"]

dumTab4 = rbind(dumTab,dumTab2)
dumTab4 = dumTab4[!is.na(bias_X1),]

dumTab4[,dumID := paste(Sim_SNPset, Sim_age, Sim_reg, sep=" - ")]
dumTab4[,dumID2 := gsub("bin","",outcome)]
dumTab4[,dumID := gsub("gamlssIA","GAM",dumID)]
dumTab4[,dumID := gsub("linMixed","LMM",dumID)]

plot5 = ggplot(dumTab4, 
               aes(x=dumID, y=bias_X1, color = dumID2, shape = outcome_type)) +
  facet_wrap(~ type,scales = "free") +
  geom_hline(yintercept = 0,color="grey") +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=bias_X1-1.96*bias_SE_X1, ymax=bias_X1+1.96*bias_SE_X1), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  #theme(axis.text.x = element_text(angle = 45)) +
  xlab("Scenario") + ylab("Bias") +
  labs(color = "Outcome", shape="Type")
plot5

filename = paste0("../results/_figures/MainFig2_Bias.png")
png(filename = filename,width = 2800, height = 1600, res=200)
print(plot5)
dev.off()


#' # Session Info ####
#' ***

sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
