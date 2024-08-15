#' ---
#' title: "Get Correlation Plots"
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
#' **Correlation Plot between raw - log - Z-scores and mean - slope - variability**
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_HPC.R")

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("2024-","24-",tag)
tag = gsub("-","",tag)

#' # Load data ####
#' ***
myFiles = list.files(path = "../results/",pattern = "02_SNPs_")

dumTab1 = foreach(i = 1:length(myFiles))%do%{
  #i=2
  loaded = load(paste0("../results/",myFiles[i]))
  myTab_loaded = get(loaded)
  ToAdd = unlist(strsplit(myFiles[i],"_"))
  ToAdd = gsub(".RData","",ToAdd)
  if(i==1){
    myTab_loaded[,setting := ToAdd[4]]
  }else{
    myTab_loaded[,setting := ToAdd[5]]
  }
  myTab_loaded
}
myTab = rbindlist(dumTab1,fill = T)

#' # Get correlation ####
#' ***
#' ## 2 Phenotypes - all settings
myPhenos = unique(myTab$model)

myTab4 = copy(myTab)[model == "combined",]
myTab4 = dcast(myTab4, SNP ~ setting, value.var=c("beta_mean","beta_slope","beta_var"))
myTab4[,beta_slope_noSlope := NULL]
myTab4[,beta_var_noVar := NULL]
myMatrix4 = as.matrix(myTab4[,-1])
CorTab4 = cor(myMatrix4,use = "pairwise.complete.obs")

myOtherPhenos = myPhenos[-1]
for(i in 1:length(myOtherPhenos)){
  #i=1
  myTab2 = copy(myTab)[model == myOtherPhenos[i]]
  
  # make long to wide
  myTab3 = dcast(myTab2, SNP ~ setting, value.var=c("beta_mean","beta_slope","beta_var"))
  
  # filter NA columns
  myTab3[,beta_slope_noSlope := NULL]
  myTab3[,beta_var_noVar := NULL]
  myMatrix = as.matrix(myTab3[,-1])
  dim(myMatrix)
  
  # correlation of everything
  CorTab = cor(myMatrix,use = "pairwise.complete.obs")
  CorTab[upper.tri(CorTab)] = CorTab4[upper.tri(CorTab4)]
  colnames(CorTab)=rep("",dim(CorTab)[1])
  rownames(CorTab)=gsub("beta_","",rownames(CorTab))
  png(file=paste0("../results/_figures/05_GenCor/CorPlot_",myOtherPhenos[i],".png"),
      width=1500,height=900,res = 200)
  corrplot(CorTab, order = "hclust",#addrect = 3,
           tl.col = "black", tl.srt = 45) 
  dev.off()

  
}  

#' ## 2 Settings - all phenotypes
#' Step 2: correlation within one model - between the phenotypes
mySettings = unique(myTab$setting)
myMainSetting = mySettings[1]
myOtherSettings = mySettings[-1]

myTab4 = copy(myTab)[setting == myMainSetting,]
myTab4 = dcast(myTab4, SNP ~ model, value.var=c("beta_mean","beta_slope","beta_var"))
myMatrix4 = as.matrix(myTab4[,-1])
CorTab4 = cor(myMatrix4,use = "pairwise.complete.obs")

for(i in 1:length(myOtherSettings)){
  #i=1
  myTab2 = copy(myTab)[setting == myOtherSettings[i]]
  
  # make long to wide
  myTab3 = dcast(myTab2, SNP ~ model, value.var=c("beta_mean","beta_slope","beta_var"))
  
  # filter NA columns
  if(i==2){
    filt = grepl("var",names(myTab3))
    myTab3 = myTab3[,!filt,with=F]
    myTab5 = copy(myTab4)[,!filt,with=F]
    myMatrix5 = as.matrix(myTab5[,-1])
    CorTab5 = cor(myMatrix5,use = "pairwise.complete.obs")
    NR_rect = 2
  }else if(i==3){
    filt = grepl("slope",names(myTab3))
    myTab3 = myTab3[,!filt,with=F]
    myTab5 = copy(myTab4)[,!filt,with=F]
    myMatrix5 = as.matrix(myTab5[,-1])
    CorTab5 = cor(myMatrix5,use = "pairwise.complete.obs")
    NR_rect = 2
  }else{
    CorTab5 = copy(CorTab4)
    NR_rect = 2
  }
  
  myMatrix = as.matrix(myTab3[,-1])
  dim(myMatrix)
  
  # correlation of everything
  CorTab = cor(myMatrix,use = "pairwise.complete.obs")
  CorTab[upper.tri(CorTab)] = CorTab5[upper.tri(CorTab5)]
  colnames(CorTab)=rep("",dim(CorTab)[1])
  rownames(CorTab)=gsub("beta_","",rownames(CorTab))
  png(file=paste0("../results/_figures/05_GenCor/CorPlot_",myOtherSettings[i],".png"),
      width=1500,height=900,res = 200)
  corrplot(CorTab, order = "hclust",#addrect = NR_rect,
           tl.col = "black", tl.srt = 45) 
  dev.off()
  
}  

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
