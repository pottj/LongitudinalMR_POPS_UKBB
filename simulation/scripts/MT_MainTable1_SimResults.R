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
#' - only the continuous outcomes
#' - only use GAMLSS
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

#' # Get table ####
#' ***
myTab2 = copy(myTab)
myTab2 = myTab2[Sim_growth=="quad"]
myTab2 = myTab2[outcome %in% c("Y2","Y3","Y4")]
myTab2 = myTab2[Sim_reg == "gamlssIA",]
myTab2[,dumID := paste(Sim_SNPset,Sim_age,sep="_")]

# I need to transform the data
myTab3 = data.table(dumID = rep(unique(myTab2$dumID),2),
                    type = rep(c("mean","slope"),each=4),
                    P_Y2 = c(myTab2[outcome=="Y2",N_sig_proc_X1],
                             myTab2[outcome=="Y2",N_sig_proc_X2]),
                    B_Y2 = c(myTab2[outcome=="Y2",bias_X1],
                             myTab2[outcome=="Y2",bias_X2]),
                    SE_BY2 = c(myTab2[outcome=="Y2",bias_SE_X1],
                               myTab2[outcome=="Y2",bias_SE_X2]),
                    P_Y3 = c(myTab2[outcome=="Y3",N_sig_proc_X1],
                             myTab2[outcome=="Y3",N_sig_proc_X2]),
                    B_Y3 = c(myTab2[outcome=="Y3",bias_X1],
                             myTab2[outcome=="Y3",bias_X2]),
                    SE_BY3 = c(myTab2[outcome=="Y3",bias_SE_X1],
                               myTab2[outcome=="Y3",bias_SE_X2]),
                    P_Y4 = c(myTab2[outcome=="Y4",N_sig_proc_X1],
                             myTab2[outcome=="Y4",N_sig_proc_X2]),
                    B_Y4 = c(myTab2[outcome=="Y4",bias_X1],
                             myTab2[outcome=="Y4",bias_X2]),
                    SE_BY4 = c(myTab2[outcome=="Y4",bias_SE_X1],
                               myTab2[outcome=="Y4",bias_SE_X2]))

myTab3[,B_Y2 := round(B_Y2,3)]
myTab3[,B_Y3 := round(B_Y3,3)]
myTab3[,B_Y4 := round(B_Y4,3)]

myTab3[,SE_BY2 := round(SE_BY2,3)]
myTab3[,SE_BY3 := round(SE_BY3,3)]
myTab3[,SE_BY4 := round(SE_BY4,3)]

myTab3[,P_Y2 := 100*P_Y2]
myTab3[,P_Y3 := 100*P_Y3]
myTab3[,P_Y4 := 100*P_Y4]

#' Change order of rows
setorder(myTab3,dumID)
setnames(myTab3,"dumID","set")

myTab3[c(1,3,5,7),set := gsub("_.*","",set)]
myTab3[c(2,4,6,8),set := gsub(".*_","",set)]

#' I will still have to make a few changes in the latex file (e.g. math mode for Y2, Y3, and Y4). 
#' 
#' # Save table ####
#' ***   
print(xtable(myTab3, type = "latex",digits = 3),file = "../results/_tables/Simulation_mainTable1.tex")

save(myTab3,file="../results/_tables/Simulation_mainTable1.RData")

excel_fn = paste0("../results/_tables/Simulation_mainTable1.xlsx")
WriteXLS("myTab3", 
         ExcelFileName=excel_fn, 
         SheetNames="Sim_main", 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

#' # Session Info ####
#' ***

sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
