#' ---
#' title: "Main simulation figure: Power heatmaps"
#' subtitle: "Main Figures"
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

outdir_results = "../results/_figures/MainFigures/"
if(dir.exists(outdir_results)==F){
  dir.create(outdir_results)
  message("Created results folder ",outdir_results, " for main figures ")
}else{
  message("Using pre-existing results folder ",outdir_results, " for main figures ")
}

#' # Get data ####
#' ***
myMainScenarios = list.files(path = "../results/_tables/",pattern = "main")
mynoSlopeScenarios = list.files(path = "../results/_tables/",pattern = "noSlope")
myScenarios = c(myMainScenarios,mynoSlopeScenarios)

dumTab = foreach(i=1:length(myScenarios))%do%{
  #i=1
  load(paste0("../results/_tables/",myScenarios[i]))
  myFlag = gsub("Simulation_complete_","",myScenarios[i])
  myFlag = gsub(".RData","",myFlag)
  myTab[,flag:=myFlag]
  myTab
}
myTab = rbindlist(dumTab, fill=T)
myTab = myTab[Sim_Y %in% c("CM2")]
myTab[,table(flag)]

#' # Get plots ####
#' ***

plotData = copy(myTab)
plotData[,dumID1 := paste(flag,Sim_X,sep="_")]
dumTab_X1 <- dcast(plotData, outcome ~ dumID1, value.var="power_X1")
dumTab_X2 <- dcast(plotData, outcome ~ dumID1, value.var="power_X2")
dumTab_X3 <- dcast(plotData, outcome ~ dumID1, value.var="power_X3")
dumTab2 = cbind(dumTab_X1,dumTab_X2[,-1],dumTab_X3[,-1])
names(dumTab2) = gsub("main_","",names(dumTab2))
names(dumTab2) = gsub("sens_","",names(dumTab2))

dumMat = as.matrix(dumTab2[,-1])
colnames(dumMat) = paste(colnames(dumMat),rep(c("X1","X2","X3"),each=12),sep=" - ")
filt1 = grepl("_X12 ",colnames(dumMat))
filt2 = grepl("_X123 ",colnames(dumMat))
filt3 = grepl("_X13 ",colnames(dumMat))

colnames(dumMat) = gsub("_X123","",colnames(dumMat))
colnames(dumMat) = gsub("_X13","",colnames(dumMat))
colnames(dumMat) = gsub("_X12","",colnames(dumMat))
rownames(dumMat) = paste0("Y",1:8)

filt = grepl("noSlope",colnames(dumMat)) & grepl("X2", colnames(dumMat))

dumMat1 = dumMat[,filt1 & !filt]
dumMat1 = dumMat1[,c(1,5,7,2,6,8,3,9,4,10)]

dumMat2 = dumMat[,filt2 & !filt]
dumMat2 = dumMat2[,c(1,5,7,2,6,8,3,9,4,10)]

dumMat3 = dumMat[,filt3 & !filt]
dumMat3 = dumMat3[,c(1,5,7,2,6,8,3,9,4,10)]

filename = paste0(outdir_results,"/Power_X12.png")
png(filename = filename,width = 2600, height = 1400, res=200)
corrplot(dumMat1, is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'), 
         col= colorRampPalette(c("#FFF5F0","#FB6A4A","#67000D"))(10),
         addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
dev.off()

filename = paste0(outdir_results,"/Power_X123.png")
png(filename = filename,width = 2600, height = 1400, res=200)
corrplot(dumMat2, is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'), 
         col= colorRampPalette(c("#FFF5F0","#FB6A4A","#67000D"))(10),
         addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
dev.off()

filename = paste0(outdir_results,"/Power_X13.png")
png(filename = filename,width = 2600, height = 1400, res=200)
corrplot(dumMat3, is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'), 
         col= colorRampPalette(c("#FFF5F0","#FB6A4A","#67000D"))(10),
         addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
dev.off()

#' # Session Info ####
#' ***

sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
