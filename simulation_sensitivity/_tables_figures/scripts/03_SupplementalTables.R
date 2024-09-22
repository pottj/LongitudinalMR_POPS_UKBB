#' ---
#' title: "Supplemental Tables for simulation"
#' subtitle: "Longitudinal MVMR"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output: html_document
#' ---
#'
#' # Introduction ####
#' ***
#' Here, I want to create the supplemental tables for the manuscript regarding the simulation study. I want in this tables
#' 
#' - one table for main simulation 
#' - one table including all sensitivity checks
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../../SourceFile_HPC.R")

#' # Get content table (tab0) ####
#' ***
{
  tab0 = data.table(Table = paste0("S",c(1:2)),
                    Title = c("Simulation results - main",
                              "Simulation results - senstivity"),
                    Source = c("simulation_main/results/_tables/",
                               "simulation_sensitivity/XXX/results/_tables/"))
  
  tab0
  
}

#' # Get Sup Tab 1 ####
#' ***
#' Main simulation
{
  loaded = load("../results/_tables/Simulation_complete_main_v1.RData")
  tab1 = get(loaded)
  names(tab1) = gsub("X1","mean",names(tab1))
  names(tab1) = gsub("X2","slope",names(tab1))
  names(tab1) = gsub("X3","var",names(tab1))
  
}

#' # Get Sup Tab 2 ####
#' ***
#' Sensitivity checks
#' 
{
  myToDoList = list.files(path = "../results/_tables/",
                          pattern = "Simulation_complete_*")
  myToDoList = myToDoList[!grepl("main",myToDoList)]
  dumTab2 = foreach(i = 1:length(myToDoList))%do%{
    #i=1
    loaded = load(paste0("../results/_tables/",myToDoList[i]))
    tab2 = get(loaded)
    names(tab2) = gsub("X1","mean",names(tab2))
    names(tab2) = gsub("X2","slope",names(tab2))
    names(tab2) = gsub("X3","var",names(tab2))
    tab2[,check := gsub("/.*","",myToDoList[i])]
    tab2
  }
  tab2 = rbindlist(dumTab2,fill=T)
  dumtab1 = copy(tab1)
  dumtab1[,check := "main"]
  tab2 = rbind(dumtab1,tab2)
  setorder(tab2,"Sim_NR")
  
}

#' # Save ####
#' ***
excel_fn = paste0("../results/SupTabs.xlsx")
WriteXLS(c("tab0","tab1","tab2"),
         ExcelFileName=excel_fn,
         SheetNames=c("Content","Tab_S1","Tab_S2"),
         AutoFilter=T,
         BoldHeaderRow=T,
         FreezeRow=1)

#' # Session Info ####
#' ***

sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
