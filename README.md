# MVMR of longitudinal exposure data

last updated: 19/03/25

[![DOI](https://zenodo.org/badge/662550343.svg)](https://doi.org/10.5281/zenodo.15052247)

This repository includes code relevant for my project **MVMR using longitudinal exposure data**. It consists of simulation studies and real data applications in POPS (Pregnancy Outcome Prediction Study)  and UK Biobank (UKB). 

## Overview

### General aim 

The aim of our study is to investigate the performance of multivariate Mendelian Randomization (MVMR) using longitudinal exposures, testing both **mean** exposure level, exposure trajectory (**slope**), and within-individual **variability** for causal effects on an outcome of interest.  

### Hypothesis 1

MVMR can separate the causal effect of mean, slope, and variability of an exposure X on an outcome Y. This will be tested in the simulation study.

### Hypothesis 2

Estimated fetal weight (EFW) has a causal effect on birth weight (BW). While this sound obvious at first, the underlying research question is if the effect is only origin from the mean EFW, the trajectory of the EFW, or also influenced by the variability of each individual. This will be tested in the POPS data, and serve as positive control. 

- Exposure: log-transformed EFW
- Outcome: 
    - 1-sample: BW in POPS (adjusted for gestational age (GA))
    - 2-sample: BW in UKB (obtained from [Neale lab](https://www.nealelab.is/uk-biobank), not adjusted for GA)

### Hypothesis 3

Total cholesterol (TC) has a causal effect on coronary artery disease (CAD). Similar as *Hypothesis 2*, this is a positive control, as we know that there is a causal effect. However, we test if the change of TC over time and the within-individual variability have an additional effect on the risk for CAD.  

- Exposure: TC
- Outcome: 
    - 1-sample: CAD in UKB (obtained from [FinnGen + panUKB meta-analysis](https://public-metaresults-fg-ukbb.finngen.fi/) "Coronary atherosclerosis" / "I9_CORATHER", though only the UKB results were used. Neale lab using linear regression for binary outcomes and is hence not considered for 1-sample MR)
    - 2-sample: CAD in [Aragam et al.](https://pubmed.ncbi.nlm.nih.gov/36474045/) (2022, latest meta-GWAS, downloaded from [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90132314))

## Structure of github repository

### helperfunctions

Various helperfunctions which I source in (**check which one of them are really necessary!**). I will try and document them as if they were within an R package (see YAML header in function for documentation).

### paper

Supplemental material for the manuscript as submitted to _Statistics in Medicine_. 

### realdata_EGG 

Scripts to run POPS analysis

1) Pre-processing 
  - SNP selection from [EGG Consortium](http://egg-consortium.org/)("Birth Weight Summary Data - Fetal GWAS (2016)") 
  - Check SNPs in POPS (AF, LD, Allele coding, ...)
  - Check genetic principal components in POPS
  - Check POPS exposure and outcome data
  
2) Get SNP associations on EWF
  - MAIN: GAMLSS with SNP-time interaction in all samples with at least 2 measurements 
  - SENS - GBR3: GAMLSS with SNP-time interaction in British ancestry samples with 3 measurements 
  - SENS - noVar: GAMLSS with SNP-time interaction in all samples with at least 2 measurements (no variability)
  - SENS - noSlope: GAMLSS without SNP-time interaction in all samples with at least 2 measurements (no slope) 

3) Get SNP associations on outcome
  - MAIN: linear regression in all samples
  - SENS: linear regression in British ancestry samples
  
4) MVMRs using all relevant combinations

5) Scripts to create figures and tables

### realdata_GLGC 

Scripts to run UKB analysis

1) Pre-processing 
  - SNP selection from [GLGC consortium data](https://csg.sph.umich.edu/willer/public/glgc-lipids2021/)("TC_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz")
  - Check SNPs in UKB (AF, LD, Allele coding, ...)
  - Check UKB exposure and outcome data
  
2) Get SNP associations on TC
  - MAIN: GAMLSS with SNP-time interaction in all samples 
  - SENS - sampleSet: GAMLSS with SNP-time interaction in samples with 
      - no statin treatment at any timepoint
      - age between 40-70
      - first measurement per year
  - SENS - noVar: GAMLSS with SNP-time interaction in all samples (no variability)
  - SENS - noSlope: GAMLSS without SNP-time interaction (no slope) 

3) MVMRs using all relevant combinations

4) Scripts to create figures and tables

### simulation

Code for all tested scenarios. In each scenario, we tested 12 settings regarding the exposure simulation:

- 3 ways to simulate the longitudinal exposure: 
    - $X_{12}$: SNPs affecting the mean and slope
    - $X_{13}$: SNPs affecting the mean and variability
    - $X_{123}$: SNPs affecting the mean, slope, and variability
- 4 causal models: causal effect of mean on outcome, $\theta_1$, either 0.3, 1.2, -1.2, or -0.3. The other causal effects are fixed with $\theta_2 = 0.3$ and $\theta_3 = 1$

## Abbreviations

- AF, allele frequency
- BW, birth weight
- CAD, coronary artery disease
- EFW, estimated fetal weight
- EGG, Early Growth Genetics Consortium
- GA, gestation age
- GAMLSS, generalized additive model for location, scale and shape
- GLGC, Global Lipids Genetics Consortium
- GWAS, genome-wide association study
- LD, linkage disequilibrium
- MR, Mendelian Randomization
- MVMR, multivariate Mendelian Randomization
- POPS, Pregnancy Outcome Prediction Study
- SNP, single nucleotide polymorphism
- TC, total cholesterol 
- UKB, UK Biobank

