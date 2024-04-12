#' ---
#' title: "Scatter Plots"
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
#' I want for each of the main EFW phenotypes scatter plots, each with four facets for mean & slope on BW and eCS. 
#' 
#' I use the suggestive significant threshold for this. (Makes plotting simpler.)
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../../SourceFile_HPC.R")

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("2024-","24-",tag)
tag = gsub("-","",tag)

#' # Get data ####
#' ***
load("../temp/04_MVMRInput_MAIN_240322.RData")
load("../results/04_MVMR_01_MAIN_240322.RData")
load("../results/01_Prep_02_LD_240317.RData")

#' # Filter data ####
#' ***
#' Filter for relevant phenotypes
myAssocs_X_long = myAssocs_X_long[grepl("efwcomb",phenotype)]
myAssocs_Y = myAssocs_Y[phenotype %in% c("pn_bw","pn_emcsall")]

MVMR_results = MVMR_results[grepl("efwcomb",exposure),]
MVMR_results = MVMR_results[outcome %in% c("pn_bw","pn_emcsall"),]
MVMR_results = MVMR_results[threshold=="suggestive",]

#' Filter for suggestive significant SNPs
myExposures = unique(myAssocs_X_long$phenotype)
myAssocs_X_long[,flag:=F]

for(i in 1:length(myExposures)){
  #i=1
  mySNPs = unique(myAssocs_X_long[pval<1e-6 & phenotype==myExposures[i],SNP])
  myAssocs_X_long[SNP %in% mySNPs & phenotype==myExposures[i],flag:=T]
  
}
myAssocs_X_long = myAssocs_X_long[flag==T,]
myAssocs_Y = myAssocs_Y[SNP %in% myAssocs_X_long$SNP,]

#' # Loop per phenotype ####
#' ***
myExposures = myAssocs_X_long[,.N,phenotype]
myExposures = myExposures[N>6,]
myExposures

dumTab1 = foreach(i=1:dim(myExposures)[1])%do%{
  #i=1
  message("Working on exposure ",myExposures[i,phenotype])
  
  # filter data
  data_GX = copy(myAssocs_X_long)
  data_GX = data_GX[phenotype==myExposures[i,phenotype],]
  data_GY = copy(myAssocs_Y)
  data_GY = data_GY[SNP %in% data_GX$SNP]
  
  # LD pruning as in MVMR function
  corTab = copy(LDTab)
  mySNPs = unique(data_GX$SNP)
  corTab= corTab[SNP1 %in% mySNPs & SNP2 %in% mySNPs,]
  corTab= corTab[value>0.1]
  corTab[,SNP2 := as.character(SNP2)]
  data_GX2 = copy(data_GX)
  setorder(data_GX2,pval)
  data_GX3 = data_GX2[!duplicated(SNP),]
  
  data_GX3[,indep := NA]
  while(sum(is.na(data_GX3$indep))!=0){
    mySNPs2 = data_GX3[is.na(indep),SNP]
    mySNP2 = mySNPs2[1]
    cor2 = corTab[SNP1 %in% mySNP2 | SNP2 %in% mySNP2]
    if(dim(cor2)[1]==0){
      data_GX3[SNP == mySNP2,indep := T]
    }else{
      mySNPs3 = unique(c(cor2$SNP1,cor2$SNP2))
      mySNPs3 = mySNPs3[!is.element(mySNPs3,mySNP2)]
      data_GX3[SNP %in% mySNPs3,indep := F]
      data_GX3[SNP == mySNP2,indep := T]
    }
  }
  table(data_GX3$indep, data_GX3$type)
  data_GX = data_GX[SNP %in% data_GX3[indep==T,SNP],]
  data_GY = data_GY[SNP %in% data_GX$SNP]
  
  # now plot
  plotData1 = copy(data_GX)
  plotData1[,myY := rep(data_GY[phenotype=="pn_bw",beta_mean],3)]
  plotData1[,myOutcome := "pn_bw"]
  
  plotData2 = copy(data_GX)
  plotData2[,myY := rep(data_GY[phenotype!="pn_bw",beta_mean],3)]
  plotData2[,myOutcome := "pn_emcsall"]
  
  plotData = rbind(plotData1,plotData2)
  
  matched = match(plotData$SNP,data_GX3$SNP)
  plotData[,myColor := data_GX3[matched,type]]
  plotData[,dumID := paste0(myOutcome," ~ ",type)]
  
  res3 = MVMR_results[exposure == myExposures[i,phenotype]]
  res3[,dumID := paste0(outcome, " ~ ", exposure_type)]
  
  myPlot = ggplot(data = plotData, aes(x = beta, y = myY, color = as.factor(myColor))) + 
    geom_point() +
    facet_wrap(~dumID,scales = "free")+        
    geom_abline(data = res3[setting=="multivariate"], aes(slope = beta_IVW,intercept=0,linetype="multivariate")) +
    geom_abline(data = res3[setting=="univariate"], aes(slope = beta_IVW,intercept=0,linetype="univariate")) +
    scale_linetype_manual("IVW estimate",values=c("univariate"=2,"multivariate"=1)) +
    labs(x=paste0("SNP effect on exposure"), 
         y=paste0("SNP effect on outcome"),
         color="Lowest pval in ")
  
  #filename1 = paste0("../figures/08_7_MVMR_top20_",flag,"_",exposure_name,"_",outcome_name,"_",tag,".tiff")
  #tiff(filename = filename1,width = 2250, height = 1125, res=200, compression = 'lzw')
  print(myPlot)
  #dev.off()
}

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))
