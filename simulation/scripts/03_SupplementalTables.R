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
#' - Table S1a: parameters - fixed over all scenarios
#' - Table S1b: parameters - variable over all scenarios
#' - Table S2: Simulation results (power, coverage, estimate, bias)
#' - Table S3: Genetic correlation in simulation 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
suppressPackageStartupMessages(library(WriteXLS))

#' # Get content table (tab0) ####
#' ***
{
  tab0 = data.table(Table = c("S1a","S1b","S2","S3"),
                    Title = c("Simulation parameters (fixed over all scenarios)", 
                              "Simulation parameters (scenario-specific)", 
                              "Simulation results (power, coverage and bias)", 
                              "Simulation results (genetic correlation)"))
  
  tab0
  
}

#' # Get Sup Tab 1 ####
#' ***
{
  tab1a = fread("../temp/parameters_fixed.txt")
  tab1b = fread("../temp/parameters_variable.txt")
  
  dummy = tab1a[1,]
  dummy[,parameter := "n_sim"]
  dummy[,value := 500]
  tab1a = rbind(dummy, tab1a)
  
  tab1b = tab1b[,-2]
  
}

#' # Get Sup Tab 2 ####
#' ***
#' 
{
  load("../result/_tables/Simulation_complete.RData")
  head(myTab)
  
  # remove unnecessary columns: N_ (always 500), empSE_SE_ (I focus on the estimate + empSE)
  myNames1 = names(myTab)
  myNames1 = myNames1[!grepl("N_",myNames1)]
  myNames1 = myNames1[!grepl("empSE_SE",myNames1)]
  
  # reorder the columns: I want the columns grouped per exposure type
  myNames2 = c(myNames1[1:4],myNames1[grepl("_mean",myNames1)],
               myNames1[grepl("_slope",myNames1)],
               myNames1[grepl("_var",myNames1)])
  colsOut<-setdiff(colnames(myTab),myNames2)
  myTab[,get("colsOut"):=NULL]
  setcolorder(myTab,myNames2)
  
  tab2 = copy(myTab)
  setnames(tab2,"Sim_NR","Scenario_NR")
  setnames(tab2,"Sim_name","Scenario_name")
}

#' # Get Sup Tab 3 ####
#' ***
#' 
{
  load("../result/_tables/Simulation_GeneticCorrelation.RData")
  head(myTab)
  
  # remove unnecessary columns: N_ (always 500), empSE_SE_ (I focus on the estimate + empSE)
  myNames1 = c("Sim_NR","Sim_name","exposure","exposure_types",
               "power","mean_correlation", "empSE")
  
  # reorder the columns: I want the columns grouped per exposure type
  colsOut<-setdiff(colnames(myTab),myNames1)
  myTab[,get("colsOut"):=NULL]
  setcolorder(myTab,myNames1)
  
  tab3 = copy(myTab)
  setnames(tab3,"Sim_NR","Scenario_NR")
  setnames(tab3,"Sim_name","Scenario_name")
}

#' # Save ####
#' ***
save(tab0,tab1a,tab1b,tab2,tab3,file="../result/_tables/SupTabs.RData")
excel_fn = paste0("../result/_tables/SupTabs.xlsx")
WriteXLS(c("tab0","tab1a","tab1b","tab2","tab3"),
         ExcelFileName=excel_fn,
         SheetNames=c("Content","Tab_S1a","Tab_S1b","Tab_S2","Tab_S3"),
         AutoFilter=T,
         BoldHeaderRow=T,
         FreezeRow=1)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
