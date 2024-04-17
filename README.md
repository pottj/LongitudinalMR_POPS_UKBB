# MVMR of longitudinal exposure data

last updated: 17/04/24

This repository includes code relevant for both the simulation study and the real data analysis using POPS (Pregnancy Outcome Prediction Study) data. 

## Overview

### General aim 

The aim of our study is to investigate the performance of multivariate Mendelian Randomization (MVMR) using longitudinal exposures, testing both **mean** exposure level, exposure trajectory (**slope**), and within-individual **variability** for causal effects on an outcome of interest.  

### Hypothesis 1

MVMR can separate the causal effect of mean, slope, and variability of an exposure X on an outcome Y. This will be tested in the simulation study.

### Hypothesis 2

Estimated fetal weight (EFW) has a causal effect on birth weight. While this sound obvious at first, the underlying research question is if the effect is only origin from the mean EFW, the trajectory of the EFW, or also influenced by the variability of each individual. This will be tested in the POPS data, and serve as positive control. 

### Hypothesis 3

The final hypothesis is whether EFW has a causal effect on emergency Cesarean Section (eCS) or not. This will only be tested if the positive control works out. 

### POPS exposure data

**Primary exposure**: estimated fetal weight (EFW) 

- absolute values (in kg)
- log-transformed values
- Z-scores (corrected for GA)
- centiles (pnorm(Z-score))

**Secondary exposures**: (not yet analyzed) 

- Linear-growth exposures (in cm)
    - Abdominal circumference (AC)
    - Femur length (FL)
    - Head circumference (HC)
    - Biparietal diameter (BPD)
- Z-score (GA corrected score, only when GA within certain ranges)
    - Abdominal circumference (AC)
    - Femur length (FL) 
    - HC/AC 
    - AC/FL

**Outcomes** 

- Birth weight (BW): used as positive control 
    - Raw values, Z-score and centile in POPS
- emergency Cesarean Section (eCS) in POPS (584 cases, 2420 controls) 

## Structure of github repository

### data

**This data will not be tracked by github!**

**Phenotypes and raw genetic data will not be stored in the repository but on the HPC rds (check source file to get path)** 

- POPS data extract as of 26/09/2023 (excel sheet, provided by Ulla Sovia)
- Various derivatives of the extracted phenotypes after QC and preprocessing
- Genetic data (filtered and mapped to samples with phenotypes)
- Data download from the GWAS Catalog on BW (date: 22/09/2023)
- PGS data download for BW (data: 13/11/2023)
- Summary statistics for BW

### helperfunctions

Various helperfunctions which I source in (**check which one of them are really necessary!**). I will try and document them as if they were within an R package (see YAML header in function for documentation).

### realdata 

Scripts to run POPS analysis

1) Pre-processing 
  - SNP selection from PGS
  - Check SNPs in POPS (AF, LD, Allele coding, ...)
  - Check genetic principal components in POPS
  - Check POPS exposure and outcome data
  
2) Get SNP associations on EWF
  - MAIN: GAMLSS with SNP-time interaction in all samples with at least 2 measurements 
  - SENS1: GAMLSS with SNP-time interaction in British ancestry samples with 3 measurements 
  - SENS2: LMM with SNP-time interaction in all samples with at least 2 measurements (no variability)
  - SENS3: GAMLSS with SNP-time interaction and time-independent sigma-function in all samples with at least 2 measurements 
  - SENS4: GAMLSS without SNP-time interaction in all samples with at least 2 measurements (no slope)
  - SENS5: GAMLSS with SNP-time interaction in all samples with at least 2 measurements and no additional covariables

3) Get SNP associations on outcome
  - MAIN: linear/logistic regression in all samples
  - SENS: linear/logistic regression in British ancestry samples
  
4) MVMRs using all relevant combinations

5) Scripts to create figures and tables

### simulation 

Test hypothesis 1 in 16 scenarios: 

- Growth: linear or quadratic
- SNP set: one shared SNP set or two distinct SNP sets for mean and slope effects
- Time: using real age or Follow-up number to simulate data and in regression model (creating/avoiding heteroscedacisity)
- Regression model: LMM or GAMLSS

### simulation v2

Run some sensitivity checks on simulation 

- MAIN: quadratic growth, age as time parameter, using GAMLSS with SNP-time interaction
- SENS1: using linear growth in regression model
- SENS2: using Follow-up number as time parameter
- SENS3: using GAMLSS without SNP-time interaction 
- SENS4: using GAMLSS with time-independent sigma-function
- SENS5: using GAMLSS with SNP-independent sigma-function

## Abbreviations

- AC, abdominal circumference
- AF, allele frequency
- BPD, biparietal diameter
- BW, birth weight
- EFW, estimated fetal weight
- FL, femur length
- GA, gestation age
- GAMLSS, generalized additive model for location, scale and shape
- GWAS, genome-wide association study
- HC, head circumference
- HPC, high performance computing
- LD, linkage disequilibrium
- LMM, linear mixed model
- MR, Mendelian Randomization
- MVMR, multivariate Mendelian Randomization
- PGS, polygenetic score
- POPS, Pregnancy Outcome Prediction Study
- QC, quality control
- rds, research data store
- SNP, single nucleotide polymorphism
- UKBB, UK Biobank

