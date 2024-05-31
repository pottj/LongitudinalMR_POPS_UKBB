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
#' Here, I want to estimate the SNP effects in a **gamlssIA** model for all relevant exposures:
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
myGD_files = myGD_files[grepl("240517",myGD_files)]
loaded1 = load(paste0(POPS_phenotypes,myGD_files[1]))
loaded1
loaded2 = load(paste0(POPS_phenotypes,myGD_files[2]))
loaded2
loaded3 = load(paste0(POPS_phenotypes,myGD_files[3]))
loaded3

myScores = list.files(path="../results/",pattern = "01_Prep_02_SNPList")
myScores = myScores[grepl("RData",myScores)]
myScore = myScores[length(myScores)]
loaded3 = load(paste0("../results/",myScore))
SNPList = SNPList_filtered

stopifnot(is.element(psam$FID,myTab_Y$POPSID))
stopifnot(pvar$ID == SNPList$SNP)
stopifnot(psam$FID == rownames(geno_mat))

pvar[,rsID := SNPList$rsid]

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
  #j=6
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

    myCols = names(data2)[c(1,2,6,32:36,39:41,43:46)]
    colsOut<-setdiff(colnames(data2),myCols)
    data2[,get("colsOut"):=NULL]
    setcolorder(data2,myCols)
    data3 = na.omit(data2)
    
    mod2 = gamlss(myX ~ myG + pn_sex + an_heightZ + an_smokstat + time +
                    (myG + pn_sex + an_heightZ + an_smokstat):time + 
                    ga2 + 
                    PC1 + PC2 + PC3 + PC4 + PC5 + random(x = as.factor(POPSID)),   
                  sigma.formula = ~myG, 
                  data = na.omit(data2), family = "NO",
                  control = gamlss.control(n.cyc = 200))
    
    dummy2 = summary(mod2)
    dummy2 = dummy2[grepl("myG",rownames(dummy2)),]
    
    res1 = data.table(regression = "gamlssIA",
                      age = "ga",
                      growth = "quad",
                      population = "all2",
                      SNP = mySNP_info$ID,
                      rsID = mySNP_info$rsID,
                      phenotype = myExposure,
                      sampleSize = length(unique(data3$POPSID)),
                      beta_mean = dummy2[1,1],
                      SE_mean = dummy2[1,2],
                      tval_mean = dummy2[1,3],
                      pval_mean = dummy2[1,4],
                      beta_slope = dummy2[2,1],
                      SE_slope = dummy2[2,2],
                      tval_slope = dummy2[2,3],
                      pval_slope = dummy2[2,4],
                      beta_var = dummy2[3,1],
                      SE_var = dummy2[3,2],
                      tval_var = dummy2[3,3],
                      pval_var = dummy2[3,4])
    res1
  }
  myAssocs_X = rbindlist(dumTab2)
  save(myAssocs_X,file=paste0("../temp/02_SNPAssocs_05_",myExposure,".RData"))
  myAssocs_X
}
myAssocs_X_gamlssIA = rbindlist(dumTab1)
myAssocs_X_gamlssIA

myAssocs_X_gamlssIA[,table(pval_mean<0.05,pval_var<0.05,phenotype)]

save(myAssocs_X_gamlssIA,file=paste0("../results/02_SNPs_05_SENS_SigmaTimeIndep.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
