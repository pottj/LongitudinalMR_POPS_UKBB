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
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(doRNG))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(gamlss))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(MendelianRandomization))

#############################
# Other helper functions 
#############################
# no helperFunctions yet ...

#############################
# path to data (only if not within the repository)
#############################
# no external data so far
