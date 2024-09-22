#############################
# this is a template source file
# please change all parameters accordingly
#############################

#############################
# R library and R packages
#############################
.libPaths()

suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(doRNG))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(gamlss))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(MendelianRandomization))
suppressPackageStartupMessages(library(pgenlibr))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(WriteXLS))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(forestploter))
suppressPackageStartupMessages(library(xtable))

#############################
# path to data (only if not within the repository)
#############################
# POPS data - restricted access only after approved project
POPS_genotyped_data = "PATH/TO/POPS/POPs_Hardcalls/"
POPS_SNP_data_EGG = "PATH/TO/POPS/POPs_Imputation/Pott/Pott_Fetal_EGG"
POPS_phenotypes = "PATH/TO/POPS_projectFiltered/POPS_phenotypes/"

# UKB genotype data - MRC BSU application data
UKB_SNP_data = "PATH/TO/UKB/genotypes-imputed/"
UKB_phenotypes = "PATH/TO/UKBB/phenotypes/ukb672224.tab"
UKB_phenotypes_filtered = "PATH/TO/UKB_projectFiltered_phenotypes/"
UKB_genotypes_filtered = "PATH/TO/UKB_projectFiltered_genotypes/"

# path to other data (GWAS summary statistics used for 2-sample approaches)
pathData = "PATH/TO/downloadedData"

