# MVMR of trajectories of fetal size on pregnancy outcomes

last updated: 24/01/24

Code relevant for preprocessing and analyzing the POPS data in my longitudinal MVMR analyses

## Overview

### General aim 

The aim of our study is to investigate the performance of different novel Mendelian Randomization (MR) techniques on longitudinal exposures, testing both mean exposure level and exposure trajectory for causal effects on an outcome of interest.  

### Hypothesis

Our hypothesis is that mean fetal size (intercept of *linMixed*, $\mu$ of *gamlss* and *gamlssIA*) and its variability during pregnancy (slope of *linMixed* and *gamlssIA*, $\sigma$ of *gamlss* and *gamlssIA*) have a causal impact on pregnancy outcomes such as delivery methods, e.g. emergency caesarean delivery (eCD). 

### Exposure & Outcomes 

**Primary exposure**: estimated fetal weight (EFW) 

- absolute values (in kg)
- Z-scores (corrected for GA)
- centiles (pnorm(Z-score))

**Secondary exposures**: 

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
    - Z-score and centile in POPS
    - Raw and inverse rank normal transformed in UKBB
- emergency caesarean delivery (eCD) in POPS (584 cases, 2420 controls) 

## Structure of github repository

### Data

#### Individual level data

**This data will not be tracked by github!**

- POPS data extract as of 26/09/2023 (excel sheet, provided by Ulla Sovia)
- Various derivatives of the extracted phenotypes after QC and preprocessing
- Genetic data (filtered and mapped to samples with phenotypes)

**Raw genetic data will not be stored here but on the HPC rds of POPS** 

#### Other data

- Data download from the GWAS Catalog on BW (date: 22/09/2023)
- PGS data download for BW (data: 13/11/2023)
- Summary statistics for BW

### Documentation

**This data will not be tracked by github!**

- Analysis plan as approved by all collaborators (pdf version) 
- Presentations as shared to all collaborators (pdf versions)

### Helperfunctions

Various helperfunctions which I source in (**check which one of them are really necessary!**). I will try and document them as if they were within an R package (see YAML header in function for documentation).

### Results

**This data is not yet tracked by github!**

Still to much chaos to do real tracking ...

Results of preprocessing are stored in data/IndividualLevelData!

### Scripts

#### 1) Preprocessing

1) Checking candidate SNP information from GWAS catalog and PGS
2) Checking POPS genotype data of candidate SNPS
3) Checking POPS PCA data
4) Checking POPS phenotype data
5) Check UKBB GWAS summary statistics for BW

#### 2) Get SNP effects

1) Get SNP effects on exposure using *linMixed* (estimating mean and slope effects)
2) Get SNP effects on exposure using *gamlss* (estimating mean and variability)
3) Get SNP effects on exposure using *gamlssIA* (estimating mean, slope and variability)
4) Get SNP effects on continuous outcome using a linear model
5) Get SNP effects on binary outcome using a generalized linear model

#### 3) Run MR analyses

1) Within-POPS (1-sample MR like), various p-value thresholds for instrument strength
2) Within-POPS (1-sample MR like), 20 best associated SNP for each exposure type (overlapping or distinct)
3) Exposure from POPS and outcome from UKBB (2-sample MR), various p-value thresholds for instrument strength
4) Exposure from POPS and outcome from UKBB (2-sample MR), 20 best associated SNP for each exposure type (overlapping or distinct)

#### 4) Evaluation

1) Shared instruments (Venn Diagrams)
2) Forest plots
3) Scatter plots
4) Tables

### Slurm

**This data will not be tracked by github!**

Some scripts needed to run via sbatch. Here I collect all used slurm calls

### Temp

**This data will not be tracked by github!**

Temporary files, which are to be deleted

## Abbreviations

- AC, abdominal circumference
- BPD, biparietal diameter
- BW, birth weight
- EFW, estimated fetal weight
- FL, femur length
- GA, gestation age
- gamlss, generalized additive model for location, scale and shape
- gamlssIA, generalized additive models for location, scale and shape with time-interaction
- GWAS, genome-wide association study
- HC, head circumference
- HPC, high performance computing
- linMixed, linear mixed model
- MR, Mendelian Randomization
- MVMR, multivariate Mendelian Randomization
- PGS, polygenetic score
- POPS, Pregnancy Outcome Prediction Study
- QC, quality control
- rds, research data store
- UKBB, UK Biobank

