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
#' - one table including all scenarios
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../../SourceFile_HPC.R")

#' # Get content table (tab0) ####
#' ***
{
  tab0 = data.table(Table = paste0("S",c(2)),
                    Title = c("Simulation results - all scenarios"),
                    Source = c("simulation/_tables_figures/results/_tables/"))
  
  tab0
  
}

#' # Get Sup Tab 2 ####
#' ***
#' 
{
  myToDoList = list.files(path = "../results/_tables/",
                          pattern = "Simulation_complete_*")
  myNames = gsub("Simulation_complete_","",myToDoList)
  myNames = gsub(".RData","",myNames)
  
  dumTab2 = foreach(i = 1:length(myToDoList))%do%{
    #i=1
    loaded = load(paste0("../results/_tables/",myToDoList[i]))
    tab2 = get(loaded)
    names(tab2) = gsub("X1","mean",names(tab2))
    names(tab2) = gsub("X2","slope",names(tab2))
    names(tab2) = gsub("X3","var",names(tab2))
    tab2[,check := myNames[i]]
    tab2
  }
  tab2 = rbindlist(dumTab2,fill=T)
  setcolorder(tab2,c(names(tab2)[56],names(tab2)[1:55]))
  setorder(tab2,"Sim_NR")
  tab2[,Sim_NR := NULL]
}

#' # Save ####
#' ***
save(tab2,file="../results/SupTabs.RData")
excel_fn = paste0("../results/SupTabs.xlsx")
WriteXLS(c("tab0","tab2"),
         ExcelFileName=excel_fn,
         SheetNames=c("Content","Tab_S2"),
         AutoFilter=T,
         BoldHeaderRow=T,
         FreezeRow=1)

#' # Session Info ####
#' ***

sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
