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
myTab_all = copy(myTab)

load("../results/_tables/Simulation_complete_F5.RData")
myTab_F5 = copy(myTab)

load("../results/_tables/Simulation_complete_F10.RData")
myTab_F10 = copy(myTab)

#' # Get figure ####
#' ***
#' ## unfiltered data 
#' Filter data 
dumTab = copy(myTab_all)
dumTab = dumTab[Sim_growth=="quad"]
dumTab = dumTab[!grepl("Y1",outcome),]

#' Create some dumID
dumTab[,dumID1 := paste(Sim_SNPset,Sim_age,Sim_reg,sep="_")]
dumTab[,dumID1 := gsub("gamlssIA","GAM",dumID1)]
dumTab[,dumID1 := gsub("linMixed","LMM",dumID1)]

dumTab_X1 <- dcast(dumTab, outcome ~ dumID1, value.var="N_sig_proc_X1")
dumTab_X2 <- dcast(dumTab, outcome ~ dumID1, value.var="N_sig_proc_X2")
dumTab_X3 <- dcast(dumTab, outcome ~ dumID1, value.var="N_sig_proc_X3")

dumTab2 = rbind(dumTab_X1,dumTab_X2)
dumTab2[,outcome := paste(rep(c("X1","X2"),each=6),outcome,sep="_")]
dumMat = as.matrix(dumTab2[,2:9])

rownames(dumMat) = gsub("_"," - ",dumTab2$outcome)
colnames(dumMat) = gsub("_"," - ",colnames(dumMat))
corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)

filename = paste0("../results/_figures/MainFig1_Power_orderedByX.png")
png(filename = filename,width = 1600, height = 1600, res=200)
corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
dev.off()

dumTab2 = rbind(dumTab_X1,dumTab_X2)
dumTab2[,outcome := paste(outcome,rep(c("X1","X2"),each=6),sep="_")]
setorder(dumTab2,outcome)
dumMat = as.matrix(dumTab2[,2:9])

rownames(dumMat) = gsub("_"," - ",dumTab2$outcome)
colnames(dumMat) = gsub("_"," - ",colnames(dumMat))
corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)

filename = paste0("../results/_figures/MainFig1_Power_orderedByY.png")
png(filename = filename,width = 1600, height = 1600, res=200)
corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
dev.off()

dumTab2 = rbind(dumTab_X1,dumTab_X2,dumTab_X3)
dumTab2[,outcome := paste(outcome,rep(c("X1","X2","X3"),each=6),sep="_")]
setorder(dumTab2,outcome)
dumMat = as.matrix(dumTab2[,2:9])

rownames(dumMat) = gsub("_"," - ",dumTab2$outcome)
colnames(dumMat) = gsub("_"," - ",colnames(dumMat))
dumMat = dumMat[!grepl("bin",rownames(dumMat)),]
corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45,na.label = " ")

filename = paste0("../results/_figures/SupFigX_Power_orderedByY_inclX3.png")
png(filename = filename,width = 1600, height = 1600, res=200)
corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45,na.label = " ")
dev.off()

#' ## filtered data (cond F-stat>5)
#' Filter data 
dumTab = copy(myTab_F5)
dumTab = dumTab[Sim_growth=="quad"]
dumTab = dumTab[!grepl("Y1",outcome),]

#' Create some dumID
dumTab[,dumID1 := paste(Sim_SNPset,Sim_age,Sim_reg,sep="_")]
dumTab[,dumID1 := gsub("gamlssIA","GAM",dumID1)]
dumTab[,dumID1 := gsub("linMixed","LMM",dumID1)]

dumTab_X1 <- dcast(dumTab, outcome ~ dumID1, value.var="N_sig_proc_X1")
dumTab_X2 <- dcast(dumTab, outcome ~ dumID1, value.var="N_sig_proc_X2")

dumTab2 = rbind(dumTab_X1,dumTab_X2)
dumTab2[,outcome := paste(outcome,rep(c("X1","X2"),each=6),sep="_")]
setorder(dumTab2,outcome)
dumMat = as.matrix(dumTab2[,2:9])
dumMat = dumMat[-c(1,2),]
dumMat = dumMat[,-c(6:8)]

rownames(dumMat) = gsub("_"," - ",dumTab2$outcome)
colnames(dumMat) = gsub("_"," - ",colnames(dumMat))
corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)

filename = paste0("../results/_figures/SupFig1_PowerF5_orderedByY.png")
png(filename = filename,width = 1600, height = 1600, res=200)
corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
dev.off()

#' ## filtered data (cond F-stat>10)
#' Filter data 
dumTab = copy(myTab_F10)
dumTab = dumTab[Sim_growth=="quad"]
dumTab = dumTab[!grepl("Y1",outcome),]

#' Create some dumID
dumTab[,dumID1 := paste(Sim_SNPset,Sim_age,Sim_reg,sep="_")]
dumTab[,dumID1 := gsub("gamlssIA","GAM",dumID1)]
dumTab[,dumID1 := gsub("linMixed","LMM",dumID1)]

dumTab_X1 <- dcast(dumTab, outcome ~ dumID1, value.var="N_sig_proc_X1")
dumTab_X2 <- dcast(dumTab, outcome ~ dumID1, value.var="N_sig_proc_X2")

dumTab2 = rbind(dumTab_X1,dumTab_X2)
dumTab2[,outcome := paste(outcome,rep(c("X1","X2"),each=6),sep="_")]
setorder(dumTab2,outcome)
dumMat = as.matrix(dumTab2[,2:9])
dumMat = dumMat[-c(1,2),]
dumMat = dumMat[,-c(6:8)]

rownames(dumMat) = gsub("_"," - ",dumTab2$outcome)
colnames(dumMat) = gsub("_"," - ",colnames(dumMat))
corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)

filename = paste0("../results/_figures/SupFig2_PowerF10_orderedByY.png")
png(filename = filename,width = 1600, height = 1600, res=200)
corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
dev.off()

#' # Session Info ####
#' ***

sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
