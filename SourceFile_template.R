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

#############################
# Other helper functions 
#############################
source("../helperfunctions/HWETest.R")

#############################
# path to data (only if not within the repository)
#############################
# POPS genotype data - restricted access only after approved project - PLINK2 format (pvar/pgen/psam)
POPS_SNP_data = "PATH/TO/POPS/IMPUTED/GENOTYPE/DATA/PROJECT_SPECIFIC_FILTERED/pfile"
