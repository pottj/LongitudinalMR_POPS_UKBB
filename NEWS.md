# NEWS file

Here I want to track my progress

## 2024/04/17

I updated the simulation: now the SNPs are more correlated (r2 = -0.6), but still not as perfect as in real data. But when I keep increasing the covariance, I cannot use mvrnorm, because the covariance matrix cannot be solved anymore. I also added several sensitivity checks (what happens when we use the wrong model?). 

The real data was not so much updated, I just added some more plots. 

I wrote the first draft and sent it out to Steve, Jessica, Ulla and Gordon (Abstract submitted to MR Conference in Bristol). I will also meet with Jessica and Marco again to discuss some open problems with gamlss (is it really better or is the default of 20 iteration just not good enough?)


## 2024/04/12

Happy easter - did not track my progress for a while ...

- included conditional F-statistics in the MVMR functions - was not pretty, because when too strong then it does not converge (not positive), and if too many instruments (real data) then it takes forever --> my solution: just use random 5% of all/nominal SNPs, then cond. F-statistics can be estimated, but they are kindof low
- I gave my BSU together talk and generate figures ...
- Tried again to get shared SNP set better simulated - currently under simulation_test (not tracked)
- refined the sensitivity tests for the real data (LMM, GBR3, sigmaTimeIndep, noTimeIA, it200, and noCovars) - unfinished because slurm script still pending

## 2024/03/22

**Simulation**

Okay, I updated the simulations and get my simulation evaluation done. But there are still some open questions

1) why is it performing so badly in shared + age for Y3 & Y4? It is kind of inverse. But I checked the code and it is alright. 
2) I thought about rerunning the simulation with checking in the gamlss analyses also mean + var, and slope + var (in addition to mean + slope + var). But how would I interpret these results? Are they helpful at all?
3) Should I run some comparison to the "classic" time-varying MVMR with get beta_GX at first time point and at the last time point and then estimate the causal effect with those 2 exposures? Would this be interesting? 

**Real data**

There are some creepy problems with the real data analysis: The gamlss algorithm does not converge, but still gives estimates. In some cases, this is obviously wrong, e.g. hc_cm in the main analysis, where almost all SNPs have pval==0. I want to run some sort of check here, to see how much I can trust the results. 

- idea 1: rerun with 200 iteration steps (default is 20) --> still in queue
- idea 2: rerun without time in the sigma model --> should then be rather similar to LMM, might be bad for the absolut values (those with increasing variance) --> still in queue

In the mean time, I just run the MVMR scripts as they are available.

## 2024/03/17

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