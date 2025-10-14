#' ---
#' title: "Get association of SNPs on TC (SENS: sample set)"
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
#' Here, I want to estimate the SNP effects in a **gamlssIA** model for TC. I will estimate estimate all the SNP effects **in a different sample set** (no statin treatment, age in 40-70, only one value per year).
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile.R")
.libPaths()

#' # Load and prep UKB data ####
#' ***
#' Load data (genetic data & phenotype data)
load(paste0(UKB_phenotypes_filtered,"/01_Prep_01_BL_FU_GP_merged_filtered.RData"))
load(paste0(UKB_genotypes_filtered,"UKB_merged.RData"))

stopifnot(is.element(myTab_long$ID,psam$FID))
stopifnot(pvar$ID == colnames(geno_mat))
stopifnot(psam$FID == rownames(geno_mat))
pvar[,rsID := ID]

#' # Get effects ####
#' ***
data1 = copy(myTab_long)
names(data1)[8:17] = gsub("_","",names(data1)[8:17])
data1 = data1[sens1==T,]

# prepare loop over SNPs
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
  
  mod2 = gamlss(TC ~ sex + myG*age + 
                  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + random(x = as.factor(ID)),   
                sigma.formula = ~ myG + sex + age, 
                data = na.omit(data2), family = "NO")
  
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

save(myAssocs_X,file=paste0("../results/02_SNPs_04_SENS_sampleSet.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
