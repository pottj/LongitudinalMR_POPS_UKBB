#' ---
#' title: "Get Supplemental Figures for real data (forest plots)"
#' subtitle: "Longitudinal MVMR in POPS"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document:
#'     toc: true
#'     number_sections: true
#'     toc_float: true
#'     code_folding: show
#' ---
#'
#' # Introduction ####
#' ***
#' **Supplemental Figures: Forest Plots**
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../../SourceFile_HPC.R")

#' # Load data ####
#' ***
load("../results/_tables/SupplementalTables_EGG.RData")

#' # Loop 1 ####
#' ***
#' I want one plot per exposure type. These plots should include 
#' 
#' - only the nominal SNP approaches
#' - 1- vs 2-sample approach
#' - uni- vs multivariable
#' - main vs sensitivity setting
#' 
myTab = copy(tab3)
myTab = myTab[threshold == "nominal_SNPs",]
myExposureTypes = c("mean","slope","variability")

myTab2 = copy(tab4)
myTab2 = myTab2[exposure_type %in% c("mean","slope_adj","var"),]
myTab2 = myTab2[threshold == "nominal_SNPs",]
setorder(myTab2,setting,exposure_type)

for(i in 1:3){
  #i=1
  plotData = copy(myTab)
  #plotData = plotData[setting == "2-sample"]
  plotData2 = copy(myTab2)
  #plotData2 = plotData2[setting == "2-sample"]
  if(i==1){
    plotData = plotData[,c(1,2,10:13)]
    myHeight = 1400
    plotData2 = plotData2[exposure_type == myExposureTypes[i]]
  }else if(i==2){
    plotData = plotData[,c(1,2,16,18:20)]
    myHeight = 1200
    plotData2 = plotData2[exposure_type == "slope_adj"]
  }else if(i==3){
    plotData = plotData[,c(1,2,22:25)]
    myHeight = 1200
    plotData2 = plotData2[exposure_type == "var"]
  }
  names(plotData)[3:6] = c("beta","SE","pval","condF")
  plotData = plotData[!is.na(beta),]
  
  plotData2 = plotData2[,c(1,2,8,9,10,11)]
  names(plotData2)[3:6] = c("beta","SE","pval","condF")
  
  plotData2[,setting2 := "MR"]
  plotData[,setting2 := "MVMR"]
  plotData = rbind(plotData,plotData2)
  
  plotData[,lowerCI95 := beta-1.96*SE]
  plotData[,upperCI95 := beta+1.96*SE]
  plotData[,Grouping := paste(setting,setting2,flag,sep=" - ")]
  setnames(plotData,"Grouping","Setting")
  plotData[,condF := round(condF,2)]
  plotData[,condF := as.character(condF)]
  setorder(plotData,flag,-setting2,-setting)
  plotData[setting!="2-sample",condF := ""]
  
  plotData$` ` <- paste(rep(" ", 50), collapse = " ")
  plotData$`Estimate \n[95% CI]` <- ifelse(is.na(plotData$SE), "",
                                         sprintf("%.2f [%.2f, %.2f]",
                                                 plotData$beta, plotData$lowerCI95, plotData$upperCI95))
  plotData[,Setting := gsub(" - MAIN","",Setting)]
  plotData[,Setting := gsub(" - SENS.*","",Setting)]
  plotData[,Setting := paste0("    ",Setting)]
  
  dummy = data.table(Setting = unique(plotData$flag))
  plotData = rbind(plotData,dummy, fill=T)
  if(i==1){
    plotData = plotData[c(21,1:4,22,5:8,23,9:12,24,13:16,25,17:20)]
  }else{
    plotData = plotData[c(17,1:4,18,5:8,19,9:12,20,13:16)]  
  }
  plotData[,Setting := gsub("SENS_","sens: ",Setting)]
  plotData[,Setting := gsub("noSlope","no slope",Setting)]
  plotData[,Setting := gsub("noVar","no variability",Setting)]
  plotData[,Setting := gsub("RIsigma","RI sigma",Setting)]
  plotData[,Setting := gsub("GBR3","data subset",Setting)]
  plotData[is.na(flag),` ` := ""]
  plotData[is.na(flag),`Estimate \n[95% CI]`:= ""]
  plotData[is.na(flag),condF := ""]
  
  dummy = plotData$Setting
  if(i==1){
    dummy[!grepl("MV",dummy)] = "#FBE3D6"
    dummy[grepl("MV",dummy)] = "#F2AA84"
  }else if(i==2){
    dummy[!grepl("MV",dummy)] = "#C2F1C8"
    dummy[grepl("MV",dummy)] = "#47D45A"
  }else if(i==3){
    dummy[!grepl("MV",dummy)] = "#CAEEFB"
    dummy[grepl("MV",dummy)] = "#61CBF4"
  }
  dummy2 = plotData$flag
  dummy[is.na(dummy2)] = "white"
  
  tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))
  
  myXlab = paste0("Causal effect of the ",myExposureTypes[i]," of EFW on BW")
  setnames(plotData,"condF","(cond) \nF-stat")
  setnames(plotData,"Setting", "Setting \n   MR approach")
  
  p2<- forest(plotData[,c(10,11,12,6)],
              est = plotData$beta,
              lower = plotData$lowerCI95, 
              upper = plotData$upperCI95,
              sizes = 0.5,
              ci_column = 2,
              ref_line = 0,
              #title = myXlab,
              theme = tm1)
  
  plot(p2)
  
  filename = paste0("../results/_figures/SupFigs/ForestPlot_",myExposureTypes[i],"_nom.png")
  png(filename = filename,width = 1700, height = myHeight, res=200)
  plot(p2)
  dev.off()
  
}

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))

