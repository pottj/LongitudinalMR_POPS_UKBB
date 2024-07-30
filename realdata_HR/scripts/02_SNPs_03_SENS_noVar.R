#' ---
#' title: "Get association of HR SNPs on exposure"
#' subtitle: "Longitudinal MVMR"
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
#' Here, I want to estimate the SNP effects in a **gamlssIA** model for HR. 
#' 
#' I will use the minimal model for the constant stage (adjust for age, sex, and phaseTime). In addition, the main model will only include a random intercept. Random slope will be added in the sensitivity checks. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_HPC.R")
.libPaths()

#' # Load and prep UKB data ####
#' ***
#' Load data (genetic data & phenotype data)
load(paste0(UKB_phenotypes_filtered,"/01_Prep_05_UKB_HR_filtered.RData"))
load(paste0(UKB_phenotypes_filtered,"/01_Prep_04_SNPData_filtered.RData"))

stopifnot(is.element(myTab_cross4$ID,psam$FID))
stopifnot(pvar2$ID == colnames(geno_mat2))
stopifnot(psam$FID == rownames(geno_mat2))
pvar2[,rsID := ID]

#' # Get effects ####
#' ***
matched = match(myTab_long4$ID,myTab_cross4$ID)
myTab_long4 = cbind(myTab_long4,myTab_cross4[matched,c(25:35,39)])
myTab_long4[Sex==0,Sex:=2]

data1 = copy(myTab_long4)
data1 = data1[stageName == "Constant" & phaseTime>30,]
badIDs = data1[speed>100,unique(ID)]
data1 = data1[!is.element(ID,badIDs)]

# prepare loop over SNPs
mySNPs = pvar2$ID
matched = match(data1$ID,psam$FID)
table(is.na(matched))

registerDoParallel(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))

dumTab2 = foreach(i = 1:length(mySNPs))%dorng%{
  #dumTab2 = foreach(i = 1:length(mySNPs))%do%{
  #i=1
  mySNP = geno_mat2[,i]
  
  data2 = copy(data1)
  data2[,myG := mySNP[matched]]
  
  mod2 = gamlss(HR ~ Sex + Age + myG*phaseTime + 
                  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + random(x = as.factor(ID)),   
                sigma.formula = ~ Sex + Age + phaseTime, 
                data = na.omit(data2), family = "NO")
  
  data3 = copy(data2)
  data3 = data3[Sex == 1,]
  mod2M = gamlss(HR ~ Age + phaseTime + myG*phaseTime + 
                   PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + random(x = as.factor(ID)),   
                 sigma.formula = ~ Age + phaseTime, 
                 data = na.omit(data3), family = "NO")
  
  data4 = copy(data2)
  data4 = data4[Sex == 2,]
  mod2F = gamlss(HR ~ Age + phaseTime + myG*phaseTime + 
                   PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + random(x = as.factor(ID)),   
                 sigma.formula = ~ Age + phaseTime, 
                 data = na.omit(data4), family = "NO")
  
  dummy2 = summary(mod2)
  dummy2 = dummy2[grepl("myG",rownames(dummy2)),]
  dummy2M = summary(mod2M)
  dummy2M = dummy2M[grepl("myG",rownames(dummy2M)),]
  dummy2F = summary(mod2F)
  dummy2F = dummy2F[grepl("myG",rownames(dummy2F)),]
  
  n1= length(unique(data2[,ID]))
  n2= length(unique(data3[,ID]))
  n3= length(unique(data4[,ID]))
  res1 = data.table(SNP = rep(pvar2[i,ID],3),
                    model = c("combined","men","women"),
                    sampleSize = c(n1,n2,n3),
                    beta_mean = c(dummy2[1,1],dummy2M[1,1],dummy2F[1,1]),
                    SE_mean =   c(dummy2[1,2],dummy2M[1,2],dummy2F[1,2]),
                    tval_mean = c(dummy2[1,3],dummy2M[1,3],dummy2F[1,3]),
                    pval_mean = c(dummy2[1,4],dummy2M[1,4],dummy2F[1,4]),
                    beta_slope = c(dummy2[2,1],dummy2M[2,1],dummy2F[2,1]),
                    SE_slope =   c(dummy2[2,2],dummy2M[2,2],dummy2F[2,2]),
                    tval_slope = c(dummy2[2,3],dummy2M[2,3],dummy2F[2,3]),
                    pval_slope = c(dummy2[2,4],dummy2M[2,4],dummy2F[2,4]))
  res1
}
myAssocs_X = rbindlist(dumTab2)
myAssocs_X

save(myAssocs_X,file=paste0("../results/02_SNPs_03_SENS_noVar.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
