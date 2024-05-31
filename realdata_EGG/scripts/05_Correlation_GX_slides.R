#' ---
#' title: "Get Correlation Plots for slides"
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
  #i=1
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

#' renaming the phenotypes
myTab[phenotype == "efwcomb", phenotype := "EFW_R"]
myTab[phenotype == "logefwcomb", phenotype := "EFW_L"]
myTab[phenotype == "efwcombZv2", phenotype := "EFW_Z"]
myTab[phenotype == "efwcombv2_cent", phenotype := "EFW_Centiles"]

#' # Get correlation ####
#' ***
myTab2 = copy(myTab)[setting == "MAIN",]
myTab2 = myTab2[phenotype != "EFW_Centiles",]

# make long to wide
myTab3 = dcast(myTab2, SNP ~ phenotype, value.var=c("beta_mean","beta_slope","beta_var"))
myMatrix = as.matrix(myTab3[,-1])
dim(myMatrix)

# correlation of everything
CorTab = cor(myMatrix,use = "pairwise.complete.obs")
colnames(CorTab)=rep("",dim(CorTab)[1])
rownames(CorTab)=gsub("beta_","",rownames(CorTab))
rownames(CorTab)=gsub("EFW_Z","Z-score",rownames(CorTab))
rownames(CorTab)=gsub("EFW_R","raw",rownames(CorTab))
rownames(CorTab)=gsub("EFW_L","log-transformed",rownames(CorTab))

rownames(CorTab)=gsub("_"," - ",rownames(CorTab))

corrplot(CorTab, 
         tl.col = "black", tl.srt = 45) 

png(file=paste0("../results/_figures/05_CorrelationPlot_Slides.png"),
    width=1500,height=1200,res = 200)
corrplot(CorTab, 
         tl.col = "black", tl.srt = 45) 
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
