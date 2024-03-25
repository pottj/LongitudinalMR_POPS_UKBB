#' ---
#' title: "Get association of PGS SNPs on outcome"
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
#' Here, I want to estimate the SNP effects in a linear model for all relevant outcomes.  
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_HPC.R")
.libPaths()

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("202.-","24-",tag)
tag = gsub("-","",tag)

#' # Load and prep PGS data ####
#' ***
myGD_files = list.files(path = POPS_phenotypes,pattern = "01_Prep_04")
loaded1 = load(paste0(POPS_phenotypes,myGD_files[1]))
loaded1
loaded2 = load(paste0(POPS_phenotypes,myGD_files[2]))
loaded2
loaded3 = load(paste0(POPS_phenotypes,myGD_files[3]))
loaded3

stopifnot(is.element(psam$FID,myTab_Y$POPSID))
stopifnot(psam$FID == rownames(geno_mat))

#' # Prepare data ####
#' ***
#' some changes in the covariables
myTab_Y[pn_sex == "MALE",pn_sex := 1]
myTab_Y[pn_sex == "FEMALE",pn_sex := 2]
myTab_Y[,pn_sex := as.numeric(pn_sex)]

myTab_Y[an_smokstat == "Never smoked",an_smokstat := 0]
myTab_Y[an_smokstat == "Quit (prepregnancy)",an_smokstat := 1]
myTab_Y[an_smokstat == "Quit (During pregnancy)",an_smokstat := 2]
myTab_Y[an_smokstat == "Currently smoking",an_smokstat := 3]
myTab_Y[,an_smokstat := as.numeric(an_smokstat)]

#' transform some of the outcomes
myTab_Y[,pn_bw := pn_bw/1000]
myTab_Y[,BW_Centile_Br1990 := BW_Centile_Br1990/100]

#' # Get effects ####
#' ***
names(myTab_Y)
myOutcomes=names(myTab_Y)[c(12,15,16,18)]
myOutcomes

dumTab1 = foreach(j=1:length(myOutcomes))%do%{
  #j=1
  myOutcome = myOutcomes[j]
  message("Working on outcome ",myOutcome)
  
  data1 = copy(myTab_Y)
  data1[,myX := get(myOutcome)]
  
  # filter data for complete data
  data1 = data1[!is.na(myX),]
  
  # prepare loop over SNPs
  mySNPs = pvar$ID
  matched = match(data1$POPSID,psam$FID)
  table(is.na(matched))
  
  registerDoParallel(4)
  
  dumTab2 = foreach(i = 1:length(mySNPs))%dopar%{
  #dumTab2 = foreach(i = 1:100)%do%{
    #i=1182
    mySNP_info = pvar[i,]
    mySNP = geno_mat[,i]
    
    data2 = copy(data1)
    data2[,myG := mySNP[matched]]
    
    myCols = c("POPSID","myX", "myG", "pn_ga_wk", "pn_sex", "an_est_age", "an_heightZ", "an_smokstat", "PC1", "PC2", "PC3", "PC4", "PC5","ancestry")
    colsOut<-setdiff(colnames(data2),myCols)
    data2[,get("colsOut"):=NULL]
    setcolorder(data2,myCols)
    data2[is.na(ancestry), ancestry := "other"]
    data3 = na.omit(data2)
    data4 = data3[ancestry == "GBR",]
    
    if(j<4){
      mod1 = lm(myX ~ myG:pn_ga_wk + pn_sex + an_est_age + an_heightZ + an_smokstat +
                  PC1 + PC2 + PC3 + PC4 + PC5, data = data3)
      mod2 = lm(myX ~ myG:pn_ga_wk + pn_sex + an_est_age + an_heightZ + an_smokstat +
                  PC1 + PC2 + PC3 + PC4 + PC5, data = data4, subset= ancestry=="GBR")
    }else{
      mod1 = glm(myX ~ myG:pn_ga_wk + pn_sex + an_est_age + an_heightZ + an_smokstat +
                 PC1 + PC2 + PC3 + PC4 + PC5, data = data3, family = "binomial")
      mod2 = glm(myX ~ myG:pn_ga_wk + pn_sex + an_est_age + an_heightZ + an_smokstat +
                 PC1 + PC2 + PC3 + PC4 + PC5, data = data4, subset= ancestry=="GBR", 
                 family = "binomial")
      
    }
    
    dummy1 = summary(mod1)$coef
    dummy1 = dummy1[grepl("myG",rownames(dummy1)),]
    
    res1 = data.table(regression = "linear",
                      age = "ga",
                      population = "all",
                      SNP = mySNP_info$ID,
                      phenotype = myOutcome,
                      sampleSize = length(unique(data3$POPSID)),
                      beta_mean = dummy1[1],
                      SE_mean = dummy1[2],
                      tval_mean = dummy1[3],
                      pval_mean = dummy1[4])
    
    
    dummy2 = summary(mod2)$coef
    dummy2 = dummy2[grepl("myG",rownames(dummy2)),]
    
    res2 = data.table(regression = "linear",
                      age = "ga",
                      population = "GBR",
                      SNP = mySNP_info$ID,
                      phenotype = myOutcome,
                      sampleSize = length(unique(data4$POPSID)),
                      beta_mean = dummy2[1],
                      SE_mean = dummy2[2],
                      tval_mean = dummy2[3],
                      pval_mean = dummy2[4])
    
    res = rbind(res1,res2)
    if(j==4) res[,regression := "logistic"]
    res
  }
  myAssocs_Y = rbindlist(dumTab2)
  myAssocs_Y
}
myAssocs_Y = rbindlist(dumTab1)
myAssocs_Y

myAssocs_Y[,table(pval_mean<0.05,phenotype,population)]

save(myAssocs_Y,file=paste0("../results/03_SNPs_02_SENS_Assocs_outcome_TIA_",tag,".RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
