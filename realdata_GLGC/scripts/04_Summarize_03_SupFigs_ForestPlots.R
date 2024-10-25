#' ---
#' title: "Get Supplemental Figures for real data (forest plots)"
#' subtitle: "Longitudinal MVMR in UKB"
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
load("../results/_tables/SupplementalTables_GLGC.RData")

#' # Loop 1 ####
#' ***
#' I want one plot per exposure type. These plots should only include the MVMR results, and only the "all SNPs setting" (there will be literally no difference between "nominal" and "all"). 
#' 
myTab = copy(tab3)
myTab = myTab[threshold == "all_SNPs",]
myExposureTypes = c("mean","slope","variability")

for(i in 1:3){
  #i=1
  plotData = copy(myTab)
  if(i==1){
    plotData = plotData[,c(1,2,10:13)]
    myHeight = 800
  }else if(i==2){
    plotData = plotData[,c(1,2,16,18:20)]
    myHeight = 700
  }else if(i==3){
    plotData = plotData[,c(1,2,22:25)]
    myHeight = 700
  }
  names(plotData)[3:6] = c("beta","SE","pval","condF")
  plotData = plotData[!is.na(beta),]
  plotData[,lowerCI95 := beta-1.96*SE]
  plotData[,upperCI95 := beta+1.96*SE]
  plotData[,Grouping := paste(setting,flag,sep=" - ")]
  plotData[,Grouping := gsub("sens_","sens: ",Grouping)]
  plotData[,Grouping := gsub("noSlope","no slope",Grouping)]
  plotData[,Grouping := gsub("noVar","no variability",Grouping)]
  plotData[,Grouping := gsub("randomEffectSigma","2 RI",Grouping)]
  setnames(plotData,"Grouping","Setting")
  plotData[,condF := round(condF,2)]
  plotData[,condF := as.character(condF)]
  plotData[setting=="2-sample",condF := ""]
  
  plotData$` ` <- paste(rep(" ", 50), collapse = " ")
  plotData$`Estimate [95% CI]` <- ifelse(is.na(plotData$SE), "",
                                      sprintf("%.2f [%.2f, %.2f]",
                                              plotData$beta, plotData$lowerCI95, plotData$upperCI95))
  
  myXlab = paste0("Causal effect of the ",myExposureTypes[i]," of TC on CAD")
  
  p2<- forest(plotData[,c(9,10,11,6)],
              est = plotData$beta,
              lower = plotData$lowerCI95, 
              upper = plotData$upperCI95,
              sizes = 0.5,
              ci_column = 2,
              ref_line = 0,
              title = myXlab)
  
  plot(p2)
  
  filename = paste0("../results/_figures/SupFigs/ForestPlot_",myExposureTypes[i],"_all_1v2sample.png")
  png(filename = filename,width = 1900, height = myHeight, res=200)
  plot(p2)
  dev.off()
  
}

#' # Loop 2 ####
#' ***
#' I want one plot per exposure type. These plots should only include the MVMR results, and only the "nominal SNPs setting" (there will be literally no difference between "nominal" and "all"). 
#' 
myTab = copy(tab3)
myTab = myTab[threshold == "nominal_SNPs",]
myExposureTypes = c("mean","slope","variability")

for(i in 1:3){
  #i=1
  plotData = copy(myTab)
  if(i==1){
    plotData = plotData[,c(1,2,10:13)]
    myHeight = 800
  }else if(i==2){
    plotData = plotData[,c(1,2,16,18:20)]
    myHeight = 700
  }else if(i==3){
    plotData = plotData[,c(1,2,22:25)]
    myHeight = 700
  }
  names(plotData)[3:6] = c("beta","SE","pval","condF")
  plotData = plotData[!is.na(beta),]
  plotData[,lowerCI95 := beta-1.96*SE]
  plotData[,upperCI95 := beta+1.96*SE]
  plotData[,Grouping := paste(setting,flag,sep=" - ")]
  plotData[,Grouping := gsub("sens_","sens: ",Grouping)]
  plotData[,Grouping := gsub("noSlope","no slope",Grouping)]
  plotData[,Grouping := gsub("noVar","no variability",Grouping)]
  plotData[,Grouping := gsub("randomEffectSigma","2 RI",Grouping)]
  setnames(plotData,"Grouping","Setting")
  plotData[,condF := round(condF,2)]
  plotData[,condF := as.character(condF)]
  plotData[setting=="2-sample",condF := ""]
  
  plotData$` ` <- paste(rep(" ", 50), collapse = " ")
  plotData$`Estimate [95% CI]` <- ifelse(is.na(plotData$SE), "",
                                         sprintf("%.2f [%.2f, %.2f]",
                                                 plotData$beta, plotData$lowerCI95, plotData$upperCI95))
  
  myXlab = paste0("Causal effect of the ",myExposureTypes[i]," of TC on CAD")
  
  p2<- forest(plotData[,c(9,10,11,6)],
              est = plotData$beta,
              lower = plotData$lowerCI95, 
              upper = plotData$upperCI95,
              sizes = 0.5,
              ci_column = 2,
              ref_line = 0,
              title = myXlab)
  
  plot(p2)
  
  filename = paste0("../results/_figures/SupFigs/ForestPlot_",myExposureTypes[i],"_nominal_1v2sample.png")
  png(filename = filename,width = 1900, height = myHeight, res=200)
  plot(p2)
  dev.off()
  
}

#' # Loop 3 ####
#' ***
#' I want one plot per exposure type. These plots should only include the MVMR results, and both "all" and the "nominal"  SNPs setting. 
#' 
myTab = copy(tab3)
myExposureTypes = c("mean","slope","variability")

for(i in 1:3){
  #i=1
  plotData = copy(myTab)
  if(i==1){
    plotData = plotData[,c(1,2,5,10:13)]
    myHeight = 1400
  }else if(i==2){
    plotData = plotData[,c(1,2,5,16,18:20)]
    myHeight = 1200
  }else if(i==3){
    plotData = plotData[,c(1,2,5,22:25)]
    myHeight = 1200
  }
  names(plotData)[4:7] = c("beta","SE","pval","condF")
  plotData = plotData[!is.na(beta),]
  plotData[,lowerCI95 := beta-1.96*SE]
  plotData[,upperCI95 := beta+1.96*SE]
  plotData[,threshold := gsub("_SNPs","",threshold)]
  plotData[,threshold := gsub("nominal","nom",threshold)]
  plotData[,threshold := gsub("all","all    ",threshold)]
  
  plotData[,Grouping := paste(setting,threshold,flag,sep=" - ")]
  plotData[,Grouping := gsub("sens_","sens: ",Grouping)]
  plotData[,Grouping := gsub("noSlope","no slope",Grouping)]
  plotData[,Grouping := gsub("noVar","no variability",Grouping)]
  plotData[,Grouping := gsub("randomEffectSigma","2 RI",Grouping)]
  setnames(plotData,"Grouping","Setting")
  plotData[,condF := round(condF,2)]
  plotData[,condF := as.character(condF)]
  plotData[setting=="2-sample",condF := ""]
  
  plotData$` ` <- paste(rep(" ", 50), collapse = " ")
  plotData$`Estimate [95% CI]` <- ifelse(is.na(plotData$SE), "",
                                         sprintf("%.2f [%.2f, %.2f]",
                                                 plotData$beta, plotData$lowerCI95, plotData$upperCI95))
  
  myXlab = paste0("Causal effect of the ",myExposureTypes[i]," of TC on CAD")
  
  p2<- forest(plotData[,c(10,11,12,7)],
              est = plotData$beta,
              lower = plotData$lowerCI95, 
              upper = plotData$upperCI95,
              sizes = 0.5,
              ci_column = 2,
              ref_line = 0,
              title = myXlab)
  
  plot(p2)
  
  filename = paste0("../results/_figures/SupFigs/ForestPlot_",myExposureTypes[i],"_1v2sample.png")
  png(filename = filename,width = 2000, height = myHeight, res=200)
  plot(p2)
  dev.off()
  
}

#' # Loop 4 ####
#' ***
#' I want a forest plot for the multi- and univariate results. I only use the 2-sample MR / MVMR results, and only the "all SNPs setting" 
#' 
myTab = copy(tab3)
myTab = myTab[threshold == "all_SNPs",]
myExposureTypes = c("mean","slope","variability")

myTab2 = copy(tab4)
myTab2 = myTab2[exposure_type %in% c("mean","slope_adj","var"),]
myTab2 = myTab2[threshold == "all_SNPs",]
setorder(myTab2,setting,exposure_type)

for(i in 1:3){
  #i=1
  plotData = copy(myTab)
  plotData = plotData[setting == "2-sample"]
  plotData2 = copy(myTab2)
  plotData2 = plotData2[setting == "2-sample"]
  if(i==1){
    plotData = plotData[,c(1,2,10:13)]
    myHeight = 800
    plotData2 = plotData2[exposure_type == myExposureTypes[i]]
  }else if(i==2){
    plotData = plotData[,c(1,2,16,18:20)]
    myHeight = 700
    plotData2 = plotData2[exposure_type == "slope_adj"]
  }else if(i==3){
    plotData = plotData[,c(1,2,22:25)]
    myHeight = 700
    plotData2 = plotData2[exposure_type == "var"]
  }
  names(plotData)[3:6] = c("beta","SE","pval","condF")
  plotData = plotData[!is.na(beta),]
  plotData2 = plotData2[,c(1,2,8,9,10,11)]
  names(plotData2)[3:6] = c("beta","SE","pval","condF")
  plotData2[,setting := "univariate"]
  plotData[,setting := "multivariate"]
  plotData = rbind(plotData,plotData2)
  
  plotData[,lowerCI95 := beta-1.96*SE]
  plotData[,upperCI95 := beta+1.96*SE]
  plotData[,Grouping := paste(setting,flag,sep=" - ")]
  plotData[,Grouping := gsub("sens_","sens: ",Grouping)]
  plotData[,Grouping := gsub("noSlope","no slope",Grouping)]
  plotData[,Grouping := gsub("noVar","no variability",Grouping)]
  plotData[,Grouping := gsub("randomEffectSigma","2 RI",Grouping)]
  setnames(plotData,"Grouping","Setting")
  plotData[,condF := round(condF,2)]
  plotData[,condF := as.character(condF)]
  plotData[setting=="2-sample",condF := ""]
  
  plotData$` ` <- paste(rep(" ", 50), collapse = " ")
  plotData$`Estimate [95% CI]` <- ifelse(is.na(plotData$SE), "",
                                         sprintf("%.2f [%.2f, %.2f]",
                                                 plotData$beta, plotData$lowerCI95, plotData$upperCI95))
  
  myXlab = paste0("Causal effect of the ",myExposureTypes[i]," of TC on CAD")
  
  p2<- forest(plotData[,c(9,10,11,6)],
              est = plotData$beta,
              lower = plotData$lowerCI95, 
              upper = plotData$upperCI95,
              sizes = 0.5,
              ci_column = 2,
              ref_line = 0,
              title = myXlab)
  
  plot(p2)
  
  filename = paste0("../results/_figures/SupFigs/ForestPlot_",myExposureTypes[i],"_all_UniMulti.png")
  png(filename = filename,width = 1900, height = myHeight, res=200)
  plot(p2)
  dev.off()
  
}

#' # Loop 5 ####
#' ***
#' I want a forest plot for the multi- and univariate results. I only use the 2-sample MR / MVMR results, and only the "nominal SNPs setting" 
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
  plotData = plotData[setting == "2-sample"]
  plotData2 = copy(myTab2)
  plotData2 = plotData2[setting == "2-sample"]
  if(i==1){
    plotData = plotData[,c(1,2,10:13)]
    myHeight = 800
    plotData2 = plotData2[exposure_type == myExposureTypes[i]]
  }else if(i==2){
    plotData = plotData[,c(1,2,16,18:20)]
    myHeight = 700
    plotData2 = plotData2[exposure_type == "slope_adj"]
  }else if(i==3){
    plotData = plotData[,c(1,2,22:25)]
    myHeight = 700
    plotData2 = plotData2[exposure_type == "var"]
  }
  names(plotData)[3:6] = c("beta","SE","pval","condF")
  plotData = plotData[!is.na(beta),]
  plotData2 = plotData2[,c(1,2,8,9,10,11)]
  names(plotData2)[3:6] = c("beta","SE","pval","condF")
  plotData2[,setting := "univariate"]
  plotData[,setting := "multivariate"]
  plotData = rbind(plotData,plotData2)
  
  plotData[,lowerCI95 := beta-1.96*SE]
  plotData[,upperCI95 := beta+1.96*SE]
  plotData[,Grouping := paste(setting,flag,sep=" - ")]
  plotData[,Grouping := gsub("sens_","sens: ",Grouping)]
  plotData[,Grouping := gsub("noSlope","no slope",Grouping)]
  plotData[,Grouping := gsub("noVar","no variability",Grouping)]
  plotData[,Grouping := gsub("randomEffectSigma","2 RI",Grouping)]
  setnames(plotData,"Grouping","Setting")
  plotData[,condF := round(condF,2)]
  plotData[,condF := as.character(condF)]
  plotData[setting=="2-sample",condF := ""]
  
  plotData$` ` <- paste(rep(" ", 50), collapse = " ")
  plotData$`Estimate [95% CI]` <- ifelse(is.na(plotData$SE), "",
                                         sprintf("%.2f [%.2f, %.2f]",
                                                 plotData$beta, plotData$lowerCI95, plotData$upperCI95))
  
  myXlab = paste0("Causal effect of the ",myExposureTypes[i]," of TC on CAD")
  
  p2<- forest(plotData[,c(9,10,11,6)],
              est = plotData$beta,
              lower = plotData$lowerCI95, 
              upper = plotData$upperCI95,
              sizes = 0.5,
              ci_column = 2,
              ref_line = 0,
              title = myXlab)
  
  plot(p2)
  
  filename = paste0("../results/_figures/SupFigs/ForestPlot_",myExposureTypes[i],"_nominal_UniMulti.png")
  png(filename = filename,width = 1900, height = myHeight, res=200)
  plot(p2)
  dev.off()
  
}

#' # Loop 6 ####
#' ***
#' I want a forest plot for the multi- and univariate results. I only use the 2-sample MR / MVMR results, and both SNP settings 
#' 
myTab = copy(tab3)
myExposureTypes = c("mean","slope","variability")

myTab2 = copy(tab4)
myTab2 = myTab2[exposure_type %in% c("mean","slope_adj","var"),]
setorder(myTab2,setting,exposure_type)

for(i in 1:3){
  #i=1
  plotData = copy(myTab)
  plotData = plotData[setting == "2-sample"]
  plotData2 = copy(myTab2)
  plotData2 = plotData2[setting == "2-sample"]
  if(i==1){
    plotData = plotData[,c(1,2,5,10:13)]
    myHeight = 1400
    plotData2 = plotData2[exposure_type == myExposureTypes[i]]
  }else if(i==2){
    plotData = plotData[,c(1,2,5,16,18:20)]
    myHeight = 1200
    plotData2 = plotData2[exposure_type == "slope_adj"]
  }else if(i==3){
    plotData = plotData[,c(1,2,5,22:25)]
    myHeight = 1200
    plotData2 = plotData2[exposure_type == "var"]
  }
  names(plotData)[4:7] = c("beta","SE","pval","condF")
  plotData = plotData[!is.na(beta),]
  plotData2 = plotData2[,c(1,2,5,8,9,10,11)]
  names(plotData2)[4:7] = c("beta","SE","pval","condF")
  
  plotData2[,setting := "univariate"]
  plotData[,setting := "multivariate"]
  plotData = rbind(plotData,plotData2)
  plotData[,threshold := gsub("_SNPs","",threshold)]
  plotData[,threshold := gsub("nominal","nom",threshold)]
  plotData[,threshold := gsub("all","all    ",threshold)]
  plotData[,setting := gsub("univariate","univariate   ",setting)]
  
  plotData[,lowerCI95 := beta-1.96*SE]
  plotData[,upperCI95 := beta+1.96*SE]
  plotData[,Grouping := paste(setting,threshold,flag,sep=" - ")]
  plotData[,Grouping := gsub("sens_","sens: ",Grouping)]
  plotData[,Grouping := gsub("noSlope","no slope",Grouping)]
  plotData[,Grouping := gsub("noVar","no variability",Grouping)]
  plotData[,Grouping := gsub("randomEffectSigma","2 RI",Grouping)]
  setnames(plotData,"Grouping","Setting")
  plotData[,condF := round(condF,2)]
  plotData[,condF := as.character(condF)]
  plotData[setting=="2-sample",condF := ""]
  
  plotData$` ` <- paste(rep(" ", 50), collapse = " ")
  plotData$`Estimate [95% CI]` <- ifelse(is.na(plotData$SE), "",
                                         sprintf("%.2f [%.2f, %.2f]",
                                                 plotData$beta, plotData$lowerCI95, plotData$upperCI95))
  
  myXlab = paste0("Causal effect of the ",myExposureTypes[i]," of TC on CAD")
  
  p2<- forest(plotData[,c(10,11,12,7)],
              est = plotData$beta,
              lower = plotData$lowerCI95, 
              upper = plotData$upperCI95,
              sizes = 0.5,
              ci_column = 2,
              ref_line = 0,
              title = myXlab)
  
  plot(p2)
  
  filename = paste0("../results/_figures/SupFigs/ForestPlot_",myExposureTypes[i],"_UniMulti.png")
  png(filename = filename,width = 2000, height = myHeight, res=200)
  plot(p2)
  dev.off()
  
}

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))

