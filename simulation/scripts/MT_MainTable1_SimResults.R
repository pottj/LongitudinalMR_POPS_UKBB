#' ---
#' title: "Evaluation POPS simulation - main table 1"
#' subtitle: "Summary of all simulations"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output: html_document
#' ---
#'
#' # Introduction ####
#' ***
#' For the paper, I want one table for the simulation results
#' 
#' - only use quadratic growth
#' - only use mean and slope
#' - only use Y2, Y3, and Y4
#' - only bias for the continuous outcomes
#' - only use GAMLSS
#' 
#' What are the columns I want? 
#' 
#' - scenario
#' - outcome
#' - bias, SE, power cont, power bin, Fstat of mean
#' - bias, SE, power cont, power bin, Fstat of slope
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
load("../results/_tables/Simulation_complete_F0.RData")

#' # Get table ####
#' ***
myTab2 = copy(myTab)
myTab2 = myTab2[Sim_growth=="quad"]
myTab2 = myTab2[outcome %in% c("Y2","Y3","Y4")]
myTab2 = myTab2[Sim_reg == "gamlssIA",]
myTab2[,dumID := paste(Sim_SNPset,Sim_age,sep="_")]

names(myTab2)
myTab2 = myTab2[,c(2,4,7, 11,12,10,10, 26,27,25,25)]
names(myTab2)[c(7,11)] = c("N_sig_proc_X1_bin","N_sig_proc_X2_bin")

myTab3 = copy(myTab)
myTab3 = myTab3[Sim_growth=="quad"]
myTab3 = myTab3[outcome %in% c("Y2bin","Y3bin","Y4bin")]
myTab3 = myTab3[Sim_reg == "gamlssIA",]

myTab2[,N_sig_proc_X1_bin := myTab3$N_sig_proc_X1]
myTab2[,N_sig_proc_X2_bin := myTab3$N_sig_proc_X2]

myTab2[,N_sig_proc_X1_bin := N_sig_proc_X1_bin * 100]
myTab2[,N_sig_proc_X2_bin := N_sig_proc_X2_bin * 100]
myTab2[,N_sig_proc_X1 := N_sig_proc_X1 * 100]
myTab2[,N_sig_proc_X2 := N_sig_proc_X2 * 100]

myTab2[,bias_X1 := round(bias_X1,3)]
myTab2[,bias_X2 := round(bias_X2,3)]

myTab2[,bias_SE_X1 := round(bias_SE_X1,3)]
myTab2[,bias_SE_X2 := round(bias_SE_X2,3)]

myTab2[,N_sig_proc_X1 := round(N_sig_proc_X1,1)]
myTab2[,N_sig_proc_X2 := round(N_sig_proc_X2,1)]
myTab2[,N_sig_proc_X1_bin := round(N_sig_proc_X1_bin,1)]
myTab2[,N_sig_proc_X2_bin := round(N_sig_proc_X2_bin,1)]

#' # Save table ####
#' ***   
print(xtable(myTab2, type = "latex",digits = 3),file = "../results/_tables/Simulation_mainTable1.tex")

save(myTab2,file="../results/_tables/Simulation_mainTable1.RData")

excel_fn = paste0("../results/_tables/Simulation_mainTable1.xlsx")
WriteXLS("myTab2", 
         ExcelFileName=excel_fn, 
         SheetNames="Sim_main", 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

#' # Session Info ####
#' ***

sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
