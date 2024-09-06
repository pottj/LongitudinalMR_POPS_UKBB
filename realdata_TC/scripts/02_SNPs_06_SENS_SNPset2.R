#' ---
#' title: "Get association of TC SNPs on exposure"
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
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_HPC.R")
.libPaths()

#' # Load and prep UKB data ####
#' ***
myList = data.table(gene = c("APOE","LDLR","CELSR2","APOC1","PCSK9",
                             "APOB","TM6SF2","CETP","HMGCR","DOCK7",
                             "TOMM40","APOC4","ABCG8","LIPC","ABO",
                             "BUD13","BACE1","APOC3","APOA4","LPL"), 
                    SNP = c("rs7412","rs73015024","rs12740374","rs141622900","rs11591147",
                            "rs6548010","rs58542926","rs183130","rs12916","rs995000",
                            "rs11668327","rs12721109","rs4245791","rs633695","rs2519093",
                            "rs964184","rs116987336","rs12718462","rs61905132","rs6993414"), 
                    varEffect = c(T,F,F,T,F,
                                  F,F,F,F,F,
                                  T,F,F,F,F,
                                  T,T,T,T,T),
                    cytoband = c("19q13.32","19p13.2" ,"01p13.3","19q13.32","01p32.3",
                                 "02p24.1" ,"19p13.11","16q13"  ,"05q13.3" ,"01p31.3",
                                 "19q13.32","19q13.32","02p21"  ,"15q21.3" ,"09q34.2",
                                 "11q23.3" ,"11q23.3", "11q23.3","11q23.3", "08p21.3"))

#' Load data (genetic data & phenotype data)
load(paste0(UKB_phenotypes_filtered,"/01_Prep_01_UKB_GP_TC.RData"))
load(paste0(UKB_phenotypes_filtered,"/01_Prep_03_SNPData.RData"))

filt = is.element(pvar$ID,myList$SNP)
table(filt)
pvar2 = pvar[filt,]
dim(geno_mat)
geno_mat2 = geno_mat[,filt]
dim(geno_mat2)
stopifnot(is.element(myTab7$ID,psam$FID))
stopifnot(pvar2$ID == colnames(geno_mat2))
stopifnot(psam$FID == rownames(geno_mat2))
pvar2[,rsID := ID]

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
  
  mod2 = gamlss(TC ~ sex + myG*exposure_age + 
                  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + random(x = as.factor(ID)),   
                sigma.formula = ~ myG + sex + exposure_age, 
                data = na.omit(data2), family = "NO")
  
  data3 = copy(data2)
  data3 = data3[sex == 1,]
  mod2M = gamlss(TC ~ myG*exposure_age +  
                   PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + random(x = as.factor(ID)),   
                 sigma.formula = ~ myG + exposure_age, 
                 data = na.omit(data3), family = "NO")
  
  data4 = copy(data2)
  data4 = data4[sex == 2,]
  mod2F = gamlss(TC ~ myG*exposure_age + 
                   PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + random(x = as.factor(ID)),   
                 sigma.formula = ~ myG + exposure_age, 
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
                    pval_slope = c(dummy2[2,4],dummy2M[2,4],dummy2F[2,4]),
                    beta_var = c(dummy2[3,1],dummy2M[3,1],dummy2F[3,1]),
                    SE_var =   c(dummy2[3,2],dummy2M[3,2],dummy2F[3,2]),
                    tval_var = c(dummy2[3,3],dummy2M[3,3],dummy2F[3,3]),
                    pval_var = c(dummy2[3,4],dummy2M[3,4],dummy2F[3,4]))
  res1
}
myAssocs_X = rbindlist(dumTab2)
myAssocs_X

save(myAssocs_X,file=paste0("../results/02_SNPs_06_SENS_SNPset2.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
