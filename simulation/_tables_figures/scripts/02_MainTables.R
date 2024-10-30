#' ---
#' title: "Main Tables for simulation"
#' subtitle: "Longitudinal MVMR"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output: html_document
#' ---
#'
#' # Introduction ####
#' ***
#' Here, I want to create the main tables for the manuscript regarding the simulation study. I want in this tables
#' 
#' - one table for each exposure $X_{12}$, $X_{13}$, and $X_{123}$
#' - all eight outcomes
#' - only $CM_2$
#' - bias, SE(bias), and power 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../../SourceFile_HPC.R")

#' # Load data ####
#' ***
load("../results/_tables/Simulation_complete_main_v1.RData")

#' # Filter data ####
#' ***
myTab = myTab[Sim_Y == "CM2",]

# get relevant columns: bias, biasSE and coverage
names(myTab)
myNames = c(names(myTab)[c(2,4)],paste0("bias_X",1:3),paste0("bias_SE_X",1:3),paste0("coverage_X",1:3))
colsOut = setdiff(colnames(myTab),myNames)
myTab[,get("colsOut"):=NULL]

# round bias and SE
myTab[,bias_X1 := round(bias_X1,3)]
myTab[,bias_X2 := round(bias_X2,3)]
myTab[,bias_X3 := round(bias_X3,3)]

myTab[,bias_SE_X1 := round(bias_SE_X1,3)]
myTab[,bias_SE_X2 := round(bias_SE_X2,3)]
myTab[,bias_SE_X3 := round(bias_SE_X3,3)]

# get power in percent
myTab[,coverage_X1 := coverage_X1*100]
myTab[,coverage_X2 := coverage_X2*100]
myTab[,coverage_X3 := coverage_X3*100]

# change column names
names(myTab) = gsub("bias","B",names(myTab))
names(myTab) = gsub("coverage","C",names(myTab))
names(myTab) = gsub("X1","mean",names(myTab))
names(myTab) = gsub("X2","slope",names(myTab))
names(myTab) = gsub("X3","var",names(myTab))

# quick check 
knitr::kable(myTab[outcome %in% c("Y2","Y5","Y6","Y8")])
knitr::kable(myTab[outcome %in% c("Y3","Y5","Y7","Y8")])
knitr::kable(myTab[outcome %in% c("Y4","Y6","Y7","Y8")])

# get tables
myTab_X12 = copy(myTab)
myTab_X12 = myTab_X12[Sim_X == "X12",]
myTab_X12[,Sim_X := NULL]

myTab_X13 = copy(myTab)
myTab_X13 = myTab_X13[Sim_X == "X13",]
myTab_X13[,Sim_X := NULL]

myTab_X123 = copy(myTab)
myTab_X123 = myTab_X123[Sim_X == "X123",]
myTab_X123[,Sim_X := NULL]

#' # Save ####
#' ***
excel_fn = paste0("../results/MainTab.xlsx")
WriteXLS(c("myTab_X12","myTab_X13","myTab_X123"),
         ExcelFileName=excel_fn,
         SheetNames=c("Sim_X12","Sim_X13","Sim_X123"),
         AutoFilter=T,
         BoldHeaderRow=T,
         FreezeRow=1)

#' # Session Info ####
#' ***

sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
