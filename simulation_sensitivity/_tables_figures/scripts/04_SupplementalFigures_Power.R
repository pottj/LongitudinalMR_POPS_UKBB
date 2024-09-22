#' ---
#' title: "Supplemental Figures for simulation - power"
#' subtitle: "Longitudinal MVMR"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output: html_document
#' ---
#'
#' # Introduction ####
#' ***
#' Here, I want to create the supplemental figures for the manuscript regarding the simulation study. I want the following scenarios:
#' 
#' - no slope
#' - no variability
#' - shared SNP set
#' - sample size as POPS
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../../SourceFile_HPC.R")

#' # Load data #### 
#' ***
myToDoList = list.files(path = "../results/_tables/",
                        pattern = "Simulation_complete_*")
myToDoList = myToDoList[!grepl("main",myToDoList)]
myToDoList = myToDoList[c(4,5,9,7)]
myToDoList

for(i in 1:length(myToDoList)){
  #i=3
  flag = gsub("Simulation_complete_","",myToDoList[i])
  flag = gsub(".RData","",flag)
  loaded = load(paste0("../results/_tables/",myToDoList[i]))
  tab2 = get(loaded)
  tab2 = tab2[Sim_Y %in% c("CM1","CM2"),]
  
  tab2[,dumID1 := paste(Sim_X,Sim_Y,sep="_")]
  
  dumTab_X1 <- dcast(tab2, outcome ~ dumID1, value.var="power_X1")
  
  if(i==1){
    dumTab_X3 <- dcast(tab2, outcome ~ dumID1, value.var="power_X3")
    dumTab2 = cbind(dumTab_X1,dumTab_X3[,-1])
    dumTab3 = copy(dumTab2)
    dumTab3 = dumTab3[,c(1,4,5,10,11,6,7,12,13,2,3,8,9),with=F]
    dumMat = as.matrix(dumTab3[,-1])
    colnames(dumMat) = paste(colnames(dumMat),rep(c("X1","X1","X3","X3"),3),sep=" - ")
    
  }else if(i==2){
    dumTab_X2 <- dcast(tab2, outcome ~ dumID1, value.var="power_X2")
    dumTab2 = cbind(dumTab_X1,dumTab_X2[,-1])
    dumTab3 = copy(dumTab2)
    dumTab3 = dumTab3[,c(1,4,5,10,11,6,7,12,13,2,3,8,9),with=F]
    dumMat = as.matrix(dumTab3[,-1])
    colnames(dumMat) = paste(colnames(dumMat),rep(c("X1","X1","X2","X2"),3),sep=" - ")
    
  }else{
    dumTab_X2 <- dcast(tab2, outcome ~ dumID1, value.var="power_X2")
    dumTab_X3 <- dcast(tab2, outcome ~ dumID1, value.var="power_X3")
    dumTab2 = cbind(dumTab_X1,dumTab_X2[,-1],dumTab_X3[,-1])
    dumTab3 = copy(dumTab2)
    dumTab3 = dumTab3[,c(1,4,5,10,11,16,17,6,7,12,13,18,19,2,3,8,9,14,15),with=F]
    dumMat = as.matrix(dumTab3[,-1])
    colnames(dumMat) = paste(colnames(dumMat),rep(c("X1","X1","X2","X2","X3","X3"),3),sep=" - ")
    
  }
  
  filename = paste0("../results/SupFig_Power_",flag,".png")
  png(filename = filename,width = 2600, height = 1400, res=200)
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'), 
           col= colorRampPalette(c("#FFF5F0","#FB6A4A","#67000D"))(10),
           addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  dev.off()
}

#' # Session Info ####
#' ***

sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")


