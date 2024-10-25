#' ---
#' title: "Get association of SNPs on TC (SENS: random effect in sigma function)"
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
#' Here, I want to estimate the SNP effects in a **gamlssIA** model for TC. 
#' 
#' Script 5: SNPs 201:250
#' 
#' Estimated time: 3:00h (SNP 1 about 80 min, first (failed) run with averaged 3.8h ... stopped after SNP 104)
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
load(paste0(UKB_phenotypes_filtered,"/01_Prep_01_UKB_GP_TC_GLGC.RData"))
load(paste0(UKB_phenotypes_filtered,"/01_Prep_04_SNPData_GLGC.RData"))

stopifnot(is.element(myTab7$ID,psam$FID))
stopifnot(pvar$ID == colnames(geno_mat))
stopifnot(psam$FID == rownames(geno_mat))
pvar[,rsID := ID]

#' # Get effects ####
#' ***
matched = match(myTab6$BSU_ID,myTab7$ID)
myTab_long = cbind(myTab6,myTab7[matched,c(2,7:16)])
myTab_long[sex==0,sex:=2]
myTab_long[,exposure_type := NULL]

data1 = copy(myTab_long)
setnames(data1,"exposure_value","TC")
setnames(data1,"BSU_ID","ID")
names(data1)[11:20] = gsub("_","",names(data1)[11:20])

# prepare loop over SNPs
pvar = pvar[201:250,]
geno_mat = geno_mat[,201:250]
mySNPs = pvar$ID
matched = match(data1$ID,psam$FID)
table(is.na(matched))

registerDoParallel(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))

dumTab2 = foreach(i = 1:length(mySNPs))%dorng%{
  #dumTab2 = foreach(i = 1:length(mySNPs))%do%{
  #i=1
  mySNP = geno_mat[,i]
  
  data2 = copy(data1)
  data2[,myG := mySNP[matched]]
  
  time1<-Sys.time()
  mod2 = gamlss(TC ~ sex + myG*exposure_age + lipLowMed +
                  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + random(x = as.factor(ID)),   
                sigma.formula = ~ myG + sex + exposure_age + lipLowMed  + random(x = as.factor(ID)), 
                data = na.omit(data2), family = "NO")
  time2<-Sys.time()
  message("\nTOTAL TIME : " ,round(difftime(time2,time1,units = "mins"),3)," minutes")
  
  dummy2 = summary(mod2)
  dummy2 = dummy2[grepl("myG",rownames(dummy2)),]
  
  n1= length(unique(data2[,ID]))
  
  res1 = data.table(SNP = rep(pvar[i,ID],1),
                    model = c("combined"),
                    sampleSize = c(n1),
                    beta_mean = c(dummy2[1,1]),
                    SE_mean =   c(dummy2[1,2]),
                    tval_mean = c(dummy2[1,3]),
                    pval_mean = c(dummy2[1,4]),
                    beta_slope = c(dummy2[2,1]),
                    SE_slope =   c(dummy2[2,2]),
                    tval_slope = c(dummy2[2,3]),
                    pval_slope = c(dummy2[2,4]),
                    beta_var = c(dummy2[3,1]),
                    SE_var =   c(dummy2[3,2]),
                    tval_var = c(dummy2[3,3]),
                    pval_var = c(dummy2[3,4]))
  res1
}
myAssocs_X = rbindlist(dumTab2)
myAssocs_X

save(myAssocs_X,file=paste0("../results/02_SNPs_05_SENS_RandomEffectSigma_SNPset5.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
