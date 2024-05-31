#' what to do? check in the real data the correlation between mean, slope and var per exposure (raw, log, Z-score, centiles)
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
myFiles = myFiles[grepl("24041",myFiles)]

dumTab1 = foreach(i = 1:length(myFiles))%do%{
  #i=2
  loaded = load(paste0("../results/",myFiles[i]))
  myTab_loaded = get(loaded)
  ToAdd = unlist(strsplit(myFiles[i],"_"))
  if(i==1){
    myTab_loaded[,setting := paste0(ToAdd[3],"_",ToAdd[4])]
  }else{
    myTab_loaded[,setting := paste0(ToAdd[3],"_",ToAdd[5])]
  }
  myTab_loaded
}
myTab = rbindlist(dumTab1,fill = T)

#' # Get correlation ####
#' ***
#' Step 1: correlation within one phenotype - between the models
myPhenos = unique(myTab$phenotype)

for(i in 1:length(myPhenos)){
  #i=1
  myTab2 = copy(myTab)[phenotype == myPhenos[i]]
  
  # make long to wide
  myTab3 = dcast(myTab2, SNP ~ setting, value.var=c("beta_mean","beta_slope","beta_var"))
  
  # filter NA columns
  myTab3[,beta_slope_05_noTimeIA := NULL]
  myTab3[,beta_var_03_LMM := NULL]
  myMatrix = as.matrix(myTab3[,-1])
  dim(myMatrix)
  
  # correlation of everything
  CorTab = cor(myMatrix,use = "pairwise.complete.obs")
  colnames(CorTab)=rep("",dim(CorTab)[1])
  rownames(CorTab)=gsub("beta_","",rownames(CorTab))
  corrplot(CorTab, type = "upper", order = "hclust", 
           tl.col = "black", tl.srt = 45,
           title =  myPhenos[i]) 
  
}  

#' Step 2: correlation within one model - between the phenotypes
mySettings = unique(myTab$setting)

for(i in 1:length(mySettings)){
  #i=1
  myTab2 = copy(myTab)[setting == mySettings[i]]
  
  # make long to wide
  myTab3 = dcast(myTab2, SNP ~ phenotype, value.var=c("beta_mean","beta_slope","beta_var"))
  
  # filter NA columns
  if(i==3){
    filt = grepl("var",names(myTab3))
    table(filt)
    myTab3 = myTab3[,!filt,with=F]
  }else if(i==5){
    filt = grepl("slope",names(myTab3))
    table(filt)
    myTab3 = myTab3[,!filt,with=F]
  }
  myMatrix = as.matrix(myTab3[,-1])
  dim(myMatrix)
  
  # correlation of everything
  CorTab = cor(myMatrix,use = "pairwise.complete.obs")
  colnames(CorTab)=rep("",dim(CorTab)[1])
  rownames(CorTab)=gsub("beta_","",rownames(CorTab))
  corrplot(CorTab, type = "upper", order = "hclust", 
           tl.col = "black", tl.srt = 45,
           title =  mySettings[i]) 
  
}  

i=1
myTab2 = copy(myTab)[setting == mySettings[i]]

# make long to wide
myTab3 = dcast(myTab2, SNP ~ phenotype, value.var=c("beta_mean","beta_slope","beta_var"))

filt = grepl("efwcombv2_cent",names(myTab3))
table(filt)
myTab3 = myTab3[,!filt,with=F]
myMatrix = as.matrix(myTab3[,-1])
dim(myMatrix)

# correlation of everything
CorTab = cor(myMatrix,use = "pairwise.complete.obs")
colnames(CorTab)=rep("",dim(CorTab)[1])
rownames(CorTab)=gsub("beta_","",rownames(CorTab))
rownames(CorTab)=gsub("efwcombZv2","EFW_Z",rownames(CorTab))
rownames(CorTab)=gsub("efwcombv2_cent","EFW_P",rownames(CorTab))
rownames(CorTab)=gsub("logefwcomb","EFW_L",rownames(CorTab))
rownames(CorTab)=gsub("efwcomb","EFW_R",rownames(CorTab))
rownames(CorTab)=gsub("_"," - ",rownames(CorTab))

corrplot(CorTab, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45) 
png(file=paste0("../results/_figures/05_CorrelationPlot.png"),
    width=1500,height=1200,res = 200)
corrplot(CorTab, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45) 
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
