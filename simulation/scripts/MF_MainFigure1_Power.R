#' ---
#' title: "Evaluation POPS simulation - main figure 1"
#' subtitle: "Summary of all simulations"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output: html_document
#' ---
#'
#' # Introduction ####
#' ***
#' For the paper, I want one figure for the power
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

#' Create some dumID
dumTab[,dumID1 := paste(Sim_SNPset,Sim_age,sep="_")]
dumTab[,dumID2 := paste(outcome,Sim_reg,sep="_")]
dumTab[,dumID2 := gsub("gamlssIA","GAM",dumID2)]
dumTab[,dumID2 := gsub("linMixed","LMM",dumID2)]

dumTab_X1 <- dcast(dumTab, dumID2 ~ dumID1, value.var="N_sig_proc_X1")
dumTab_X2 <- dcast(dumTab, dumID2 ~ dumID1, value.var="N_sig_proc_X2")

dumTab2 = cbind(dumTab_X1,dumTab_X2[,c(2:5)])
names(dumTab2)[2:9] = paste(names(dumTab2)[2:9],rep(c("X1","X2"),each=4),sep="_") 
dumMat = as.matrix(dumTab2[,2:9])

rownames(dumMat) = gsub("_"," - ",dumTab2$dumID2)
colnames(dumMat) = gsub("_"," - ",colnames(dumMat))
corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)

filename = paste0("../results/_figures/MainFig1_Power.png")
png(filename = filename,width = 1600, height = 1600, res=200)
corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
dev.off()


#' # Session Info ####
#' ***

sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
