#' ---
#' title: "Get association of PGS SNPs on exposure"
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
#' Here, I want to estimate the SNP effects in a **linear mixed** model for all relevant exposures:
#' 
#' - Population: **all2** 
#' - Age: **ga** 
#' - Growth: **quadratic** 
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
#' Load genetic data (data from the PGS, not from GWAS Catalog)
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
myTab_X[pn_sex == "MALE",pn_sex := 1]
myTab_X[pn_sex == "FEMALE",pn_sex := 2]
myTab_X[,pn_sex := as.numeric(pn_sex)]

#' transform some of the exposures
myTab_X[,efwcomb := efwcomb/1000]
myTab_X[,efwcombv2_cent := efwcombv2_cent/100]
myTab_X[,logefwcomb := log(efwcomb)]

#' add maternal height, maternal age, maternal smoking status, and maternal gestational diabetes status
matched = match(myTab_X$POPSID,myTab_Y$POPSID)
myTab_X[,an_est_age := myTab_Y[matched,an_est_age]]
myTab_X[,an_heightZ := myTab_Y[matched,an_heightZ]]
myTab_X[,an_smokstat := myTab_Y[matched,an_smokstat]]
myTab_X[,pn_gdm := myTab_Y[matched,pn_gdm_diet_medication]]

myTab_X[an_smokstat == "Never smoked",an_smokstat := 0]
myTab_X[an_smokstat == "Quit (prepregnancy)",an_smokstat := 1]
myTab_X[an_smokstat == "Quit (During pregnancy)",an_smokstat := 2]
myTab_X[an_smokstat == "Currently smoking",an_smokstat := 3]
myTab_X[,an_smokstat := as.numeric(an_smokstat)]

myTab_X[,ga2 := ga^2]

#' # Get effects ####
#' ***
names(myTab_X)
myExposures=names(myTab_X)[c(18,38,23,24)]
myExposures

dumTab1 = foreach(j=1:length(myExposures))%do%{
  #j=5
  myExposure = myExposures[j]
  message("Working on exposure ",myExposure)
  
  data1 = copy(myTab_X)
  data1[,myX := get(myExposure)]
  
  # filter data for 2 or more time points per sample
  dummy = data1[!is.na(myX),.N,by=POPSID]
  data1 = data1[POPSID %in% dummy[N>=2,POPSID]]
  #data1 = data1[POPSID %in% dummy[N==3,POPSID]]
  #data1 = data1[!is.na(ancestry),]
  
  # set time parameter
  # data1[,time := as.numeric(scan)]
  data1[,time := ga]
  # data1[,time := scale(ga)]
  # data1[,ga2 := scale(ga2)]
  
  # prepare loop over SNPs
  mySNPs = pvar$ID
  matched = match(data1$POPSID,psam$FID)
  table(is.na(matched))
  
  registerDoParallel(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))
  
  dumTab2 = foreach(i = 1:length(mySNPs))%dorng%{
    #dumTab2 = foreach(i = 1:20)%do%{
    #i=1
    mySNP_info = pvar[i,]
    mySNP = geno_mat[,i]
    
    data2 = copy(data1)
    data2[,myG := mySNP[matched]]
    
    mod1 = lmer(myX ~ myG + pn_sex + an_heightZ + an_smokstat + time +
                  (myG + pn_sex + an_heightZ + an_smokstat):time + 
                  ga2 + 
                  (1|POPSID) + PC1 + PC2 + PC3 + PC4 + PC5, data = data2)
    
    dummy1 = summary(mod1)$coef
    dummy1 = dummy1[grepl("myG",rownames(dummy1)),]
    
    res1 = data.table(regression = "linMixed",
                      age = "ga",
                      growth = "quad",
                      population = "all2",
                      SNP = mySNP_info$ID,
                      phenotype = myExposure,
                      sampleSize = length(unique(mod1@frame$POPSID)),
                      beta_mean = dummy1[1,1],
                      SE_mean = dummy1[1,2],
                      tval_mean = dummy1[1,3],
                      beta_slope = dummy1[2,1],
                      SE_slope = dummy1[2,2],
                      tval_slope = dummy1[2,3])
    res1[,pval_mean := pnorm(-abs(beta_mean/SE_mean))*2]
    res1[,pval_slope := pnorm(-abs(beta_slope/SE_slope))*2]
    res1 = res1[,c(1:10,14,11:13,15)]
    res1
  }
  myAssocs_X = rbindlist(dumTab2)
  save(myAssocs_X,file=paste0("../temp/02_SNPAssocs_03_",myExposure,".RData"))
  myAssocs_X
}
myAssocs_X_linearMixed = rbindlist(dumTab1)
myAssocs_X_linearMixed

myAssocs_X_linearMixed[,table(pval_mean<0.05,pval_slope<0.05,phenotype)]

save(myAssocs_X_linearMixed,file=paste0("../results/02_SNPs_03_SENS_LMM_",tag,".RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
