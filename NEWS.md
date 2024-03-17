# NEWS file

Here I want to track my progress

## 17/03/2024

Updated the repositories structure

Updated simulation

- only "absolute values" and scenarios for 
    - shared vs. distinct SNP sets
    - age as actual age vs. scan number
    - growth linear vs. quadratic
    - regression model gamlssIA vs. linMixed

Update real data
    - rerun preparation steps (only POPS data)
    - started sbatch for exposures

## Older notes (unsorted)



1)    DONE get the MVMR with GA adjustment - add EGG_GA data as additional exposure and use in model 
2)    DONE get the MVMR with cesarean section from Sakaue (UKB based, but "born by CS"), just to test



3) get simulation done the way I need them - new repository? part of this one?
    - I want to check the correlation, but
      in my previous simulations I create the allele scores to be independent! Hence there is no correlation
      But how plausible is this? I mean, is it realistic to consider a SNP with a positive effect on the mean, but with an negative effect on the slope? Is there no nicer way to create two independent scores? Should I do it the same way I simulate the random effects? 
      
      I am confused again ... why on earth is simulation stuff so difficult???
      
      
      what is realistic? check my SNP data from POPS and compare effect sizes for mean and slope (using SNPs with MAF within 0.1 and 0.4) --> mean(beta_mean) and mean(beta_slope) 
      
      
      
      
      
      12/02/2024
      new plan
      simulation 1: simple linear growth, simulate data for 0-40, select 3 time points (age), and then test the 3 regression models linMixed, gamlss, and gamlssIA and the 3 exposure A, Z, and C
          - fix shared SNP set
          - fix 3 time points
          
      simulation 2: add ga^2 to simple linear growth model and repeat everything