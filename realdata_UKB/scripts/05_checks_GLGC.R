#' ---
#' title: "Check MR TC on CAD"
#' subtitle: "Longitudinal MVMR in UKB"
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
#' I want to check how well the TC estimates from the gamlss correlate to the GLGC data. Also, I want to check what the univariate effect in GLGC would have been, using either all SNPs, or the 44 selected SNPs in the sensitivity run. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_HPC.R")

#' # Load data ####
#' ***
load("../results/02_SNPs_01_MAIN.RData")
load("../results/01_Prep_03_SNPList.RData")

stopifnot(myAssocs_X$SNP == SNPList$rsID)

myAssocs_X[,flag := F]
myAssocs_X[pval_mean<5e-8 & pval_var>0.05,flag := T]
myAssocs_X[pval_var<5e-8 & pval_mean>0.05,flag := T]
myAssocs_X[,table(flag)]

filt = myAssocs_X$flag

#' # Check correlation ####
#' ***
plot(myAssocs_X$beta_mean, SNPList$GLGC_beta)
abline(0,1)
cor.test(myAssocs_X$beta_mean, SNPList$GLGC_beta)

plot(myAssocs_X$beta_mean[filt], SNPList$GLGC_beta[filt])
abline(0,1)
cor.test(myAssocs_X$beta_mean[filt], SNPList$GLGC_beta[filt])

plot(myAssocs_X$beta_slope, SNPList$GLGC_beta)
abline(0,1)
cor.test(myAssocs_X$beta_slope, SNPList$GLGC_beta)

plot(myAssocs_X$beta_slope[filt], SNPList$GLGC_beta[filt])
abline(0,1)
cor.test(myAssocs_X$beta_slope[filt], SNPList$GLGC_beta[filt])

plot(myAssocs_X$beta_var, SNPList$GLGC_beta)
abline(0,1)
cor.test(myAssocs_X$beta_var, SNPList$GLGC_beta)

plot(myAssocs_X$beta_var[filt], SNPList$GLGC_beta[filt])
abline(0,1)
cor.test(myAssocs_X$beta_var[filt], SNPList$GLGC_beta[filt])

#' # Check MR of TC on CAD ####
#' ***
mr_obj = mr_input(bx = SNPList$GLGC_beta,
                  bxse = SNPList$GLGC_SE,
                  by = SNPList$Aragam_beta,
                  byse = SNPList$Aragam_SE,
                  exposure = "TC (GLGC)",
                  outcome = "CAD (Aragam)",)
mr_ivw(mr_obj)

mr_obj = mr_input(bx = SNPList$GLGC_beta,
                  bxse = SNPList$GLGC_SE,
                  by = SNPList$UKB_beta,
                  byse = SNPList$UKB_SE,
                  exposure = "TC (GLGC)",
                  outcome = "CAD (UKB)",)
mr_ivw(mr_obj)

mr_obj = mr_input(bx = SNPList$GLGC_beta[filt],
                  bxse = SNPList$GLGC_SE[filt],
                  by = SNPList$Aragam_beta[filt],
                  byse = SNPList$Aragam_SE[filt],
                  exposure = "TC (GLGC)",
                  outcome = "CAD (Aragam)",)
mr_ivw(mr_obj)

mr_obj = mr_input(bx = SNPList$GLGC_beta[filt],
                  bxse = SNPList$GLGC_SE[filt],
                  by = SNPList$UKB_beta[filt],
                  byse = SNPList$UKB_SE[filt],
                  exposure = "TC (GLGC)",
                  outcome = "CAD (UKB)",)
mr_ivw(mr_obj)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
