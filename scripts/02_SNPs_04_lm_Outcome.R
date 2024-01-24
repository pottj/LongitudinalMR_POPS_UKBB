#' ---
#' title: "Get association of PGS SNPs on outcome (3)"
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

#' transform some of the outcomes
hist(myTab$BW_SDS_Br1990)

myTab[,pn_bwpct := pn_bwpct/100]
myTab[,pn_bw_cust_v678 := pn_bw_cust_v678/100]
myTab[,BW_Centile_Br1990 := BW_Centile_Br1990/100]


#' # Get effects ####
#' ***
names(myTab)
myOutcomes=names(myTab)[c(12:16)]

dumTab1 = foreach(j=1:length(myOutcomes))%do%{
  #j=5
  myOutcome = myOutcomes[j]
  message("Working on outcome ",myOutcome)
  
  data1 = copy(myTab)
  data1[,myX := get(myOutcome)]
  
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
    
    mod1 = lm(myX ~ myG + pn_ga_wk + pn_sex + an_est_age + an_height + an_smokstat + 
                  PC_01 + PC_02 + PC_03 + PC_04 + PC_05, data = data2)
    
    dummy1 = summary(mod1)$coef
    dummy1 = dummy1[grepl("myG",rownames(dummy1)),]
    
    res1 = data.table(model = "linear",
                      SNP = mySNP_info$ID,
                      phenotype = myOutcome,
                      beta_mean = dummy1[1],
                      SE_mean = dummy1[2],
                      tval_mean = dummy1[3],
                      pval_mean = dummy1[4])
    res1
  }
  myAsscos_X = rbindlist(dumTab2)
  #save(myAsscos_X,file=paste0("../temp/08_3_Associations_SNPs_outcome_lin_",myOutcome,"_",tag,".RData"))
  myAsscos_X
}
myAsscos_X_linear = rbindlist(dumTab1)
myAsscos_X_linear

myAsscos_X_linear[model=="linear",table(pval_mean<0.05,phenotype)]

save(myAsscos_X_linear,file=paste0("../results/08_3_Associations_SNPs_outcome_lin_",tag,".RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
