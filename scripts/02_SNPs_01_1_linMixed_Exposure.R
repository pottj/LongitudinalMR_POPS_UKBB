#' ---
#' title: "Get association of PGS SNPs on exposure (1)"
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
#' Here, I want to estimate the SNP effects in a linear mixed model for all relevant exposures. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_HPC.R")
.libPaths()

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("202.-","23-",tag)
tag = gsub("-","",tag)

#' # Load and prep PGS data ####
#' ***
load("../data/IndividualLevelData/06_Outcome_231126.RData")
myTab

load("../data/IndividualLevelData/03_LongitudinalExposure_filtered_231110.RData")
myTab_X

pvar1 = NewPvar("../../rds/rds-obsgynae-POPs-DC0Q0rbedlk/POPs_Imputation/Pott/Pott_Fetal_PGS.pvar")
pgen = NewPgen("../../rds/rds-obsgynae-POPs-DC0Q0rbedlk/POPs_Imputation/Pott/Pott_Fetal_PGS.pgen", pvar=pvar1)
pvar = fread("../../rds/rds-obsgynae-POPs-DC0Q0rbedlk/POPs_Imputation/Pott/Pott_Fetal_PGS.pvar")
psam = fread("../../rds/rds-obsgynae-POPs-DC0Q0rbedlk/POPs_Imputation/Pott/Pott_Fetal_PGS.psam")

geno_mat <- ReadList(pgen, variant_subset = c(1:dim(pvar)[1]) , meanimpute=F)
dim(geno_mat)

#' # Prepare data ####
#' ***
#' Match genetic data
matched = match(myTab_X$POPSID, myTab$POPSID)
table(is.na(matched))
myTab2 = copy(myTab)
myTab2 = myTab2[matched,]
myTab_X = cbind(myTab_X,myTab2[,c(2,3,8,25:36),with=F])

#' some changes in the covariables
myTab_X[pn_sex == "MALE",pn_sex := 1]
myTab_X[pn_sex == "FEMALE",pn_sex := 2]
myTab_X[,pn_sex := as.numeric(pn_sex)]

#' Reformate genotype matrix 
dim(geno_mat)
rownames(geno_mat)[1:10]
colnames(geno_mat)[1:10]
rownames(geno_mat) = psam$`#FID`
colnames(geno_mat) = pvar$ID
class(geno_mat)

table(is.element(psam$`#FID`,myTab_X$POPSID))
filt = is.element(psam$`#FID`,myTab_X$POPSID)
geno_mat = geno_mat[filt,]
dim(geno_mat)
psam = psam[filt,]

table(psam$`#FID` == myTab_X$POPSID)
matched = match(psam$`#FID`,myTab_X$POPSID  )
table(is.na(matched))
table(psam$SEX,myTab_X[matched,pn_sex])

#' remove one outlier
myTab_X[efw_zscoreHadlock>5,]
myTab_X = myTab_X[POPSID != 211,]
myTab = myTab[POPSID != 211]

#' transform some of the exposures
myTab_X[,efw := efw/1000]
myTab_X[,efw2 := efw2/1000]
myTab_X[,efwcomb := efwcomb/1000]

myTab_X[,efw_centileHadlock := efw_centileHadlock/100]
myTab_X[,efw_centileHadlock2 := efw_centileHadlock2/100]
myTab_X[,efwcombv2_cent := efwcombv2_cent/100]


#' # Get effects ####
#' ***
names(myTab_X)
myExposures=names(myTab_X)[c(12:26,29,30)]

dumTab1 = foreach(j=1:length(myExposures))%do%{
  #j=5
  myExposure = myExposures[j]
  message("Working on exposure ",myExposure)
  
  data1 = copy(myTab_X)
  data1[,myX := get(myExposure)]
  
  # create some plot
  # ggp1  = ggplot(data1, aes(x=ga, y=myX, col=as.factor(pn_sex))) +
  #   geom_point()+
  #   labs(x="Gestational Week",y=myExposure, color="Babys sex") +
  #   theme(legend.position = "none") + theme_classic() +
  #   geom_smooth(method = "loess",
  #               formula = y ~ x)
  # print(ggp1)
  
  # prepare loop over SNPs
  mySNPs = pvar$ID
  matched = match(data1$POPSID,psam$`#FID`)
  table(is.na(matched))
  
  registerDoParallel(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))

  dumTab2 = foreach(i = 1:length(mySNPs))%dorng%{
  #dumTab2 = foreach(i = 1:100)%do%{
      #i=1
    mySNP_info = pvar[i,]
    mySNP = geno_mat[,i]

    data2 = copy(data1)
    data2[,myG := mySNP[matched]]

    mod1 = lmer(myX ~ myG + ga + (1|POPSID) + myG:ga + pn_sex + an_est_age + an_height + an_smokstat + 
                  PC_01 + PC_02 + PC_03 + PC_04 + PC_05, data = data2)

    dummy1 = summary(mod1)$coef
    dummy1 = dummy1[grepl("myG",rownames(dummy1)),]
    
    res1 = data.table(model = "linearMixed",
                      SNP = mySNP_info$ID,
                      phenotype = myExposure,
                      beta_mean = dummy1[1,1],
                      SE_mean = dummy1[1,2],
                      tval_mean = dummy1[1,3],
                      beta_var = dummy1[2,1],
                      SE_var = dummy1[2,2],
                      tval_var = dummy1[2,3])
    res1[,pval_mean := pnorm(-abs(beta_mean/SE_mean))*2]
    res1[,pval_var := pnorm(-abs(beta_var/SE_var))*2]
    res1 = res1[,c(1:6,10,7:9,11)]
    res1
  }
  myAsscos_X = rbindlist(dumTab2)
  #save(myAsscos_X,file=paste0("../temp/08_1_Associations_SNPs_exposure_linMix_",myExposure,"_",tag,".RData"))
  myAsscos_X
}
myAsscos_X_linearMixed = rbindlist(dumTab1)
myAsscos_X_linearMixed

myAsscos_X_linearMixed[model=="linearMixed",table(pval_mean<0.05,pval_var<0.05,phenotype)]

save(myAsscos_X_linearMixed,file=paste0("../results/08_1_Associations_SNPs_exposure_linMix_",tag,".RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
