# NEWS file

Here I want to track my progress

## 2025/03/19

- submit manuscript to _Statistics in Medicine_

## 2024/11/05

- Found an error in the sens_A_GxE scripts

## 2024/10/25

- run analysis with random intercept in sigma model 
    - EGG: no difference
    - UKBB: some difficulties because of run time - split the SNP set into 7 subsets and run in parallel was easier than getting more cores / nodes / time when using 1 scripts. But at the end it worked. No big difference (just when using "all" SNPs, and not the signigicant ones)
    - simulation: chaos, because it seems to take forever. I reduced to 10 simulations each. Still have a run time of 15 hours. Some weird problems with scenario 5, which doesn't make sense as the code is the same as for scenario 6, which finished after 13 h (after 36 h scenario 5 still not finished, just reached SNP 3) - workaround in place, but I want scenario 5 to run through as well. Maybe increase to 50 simulations? Should not affect run time as I could increase the requested cores? (I don't expect a big change here, but 50 sounds nicer than 10 simulations)
    
## 2024/10/11

- welcome back from holiday
- updating simulation: run some scenarios with 
    - fixed genetic correlation between mean-variability and slope-variability
    - gene-environment interaction (check pleiotropy)
    
## 2024/09/20

- updated real data analyses: after feedback from Steve, I streamlined the analysis plan for the real data, to keep it focused:
    - strict SNP selection procedure (EGG for EFW; GLGC for TC)
    - I only test one exposure each (EFW log-transformed and TC)
    - I only test two outcomes each (1-sample: BW or CAD within POPS or UKB, respectively; 2-sample: BW or CAD using UKB (Neale lab) or Aragam et al. summary statistics)
    - sensitivity tests: no variability, no slope, different sample set. For TC only: different SNP set (in POPS I did not have these kind of SNP sets)
- archived old versions: 
    - using PGS in POPS (not enough SNPs, better to use some consortia data)
    - using HR in UKB (longitudinal in seconds, not meaningful)
    - using TC in UKB (SNP selection based on trajGWAS results, probably biased)
- moved some other stuff into test area: 
    - simulation: with random intercept and random slope -> still problems with positive definit matrix in GAMLSS
    - real data: TC sex-stratified (and women also by menopause) -> follow-up project

## 2024/09/13

- check CAD in UKB: there were really small effect sizes, so I rerun it in our BSU data again (including sex-stratified). Lower power than the data from panUKB, but high correlation between the estimates (0.9)
- plotting 
- rerunning some stuff in POPS for paper

## 2024/09/06

- finish main part of TC analyses
- continue paper writing

## 2024/08/21

- added TC analysis
- finished simulation stuff
- paper writing
- holiday :)

## 2024/08/02

Long time testing the random slope model, but it still only works for some seeds. Problem is that for some SNPs the covariance matrix is not invertible, e.g. one covariable is linear combination of another. Not sure why this is only a problem when including random slope. 

I also updated the main simulation to create stronger instruments here. I simply changed the factor for the allele score from 1 to 0.5. The simulations are still running. 

I updated the evaluation script in extract coverage and empSE as well. 

The runs on the HR data is finished, but I want to recheck the phase selection: in the constant phase, I only want data points after 30 seconds to be sure that we have "heart rate in movement". Also, for the ramp-up, I simply want protocol ID 66 for women and 88 for men. In both cases, the ramp-up is a factor 3, from 30 or 40 to 90 or 120 watt.


## 2024/07/05

I added the data prep for UKB HR data and did the first non-genetic model tests. I spoke with Marco, might be best to try the fully adjusted model including sex-IA, as the protocol might create such sex-interaction with the covariables. 

Next week: run GAMLSS and download matching CAD summary statistics for the MVMR. 

## 2024/06/28

Two weeks in one, because I was at the MR Conference in Bristol, and had some days off. 

I did some tests regarding longitudinal data from the UKBB (heart rate during bicycle exercise), but these scripts are not yet staged (still in test folder). This should be added next week. 

So far, simulation and POPS analyses are complete. Back to writing. 

## 2024/06/14

Okay, what did I do in the last week? 

### Simulation

- Run sensitivity simulation with binary outcomes with prevalence of 0.25 (instead of 0.5, more similar to POPS)
- Rerun simulation evaluation as there was a typo in my code (old $/theta$s)

### Real data

- add a script for a nicer forest plot for the 1-sample and 2-sample MVMR results

## 2024/06/07

Okay, what did I do in the last week? 

### Simulation

**Feedback Steve**: maybe focus on different values for $\theta_1$? 

- CM1: $\theta_1 = 0.3 = \theta_2$ (as before)
- CM2: $\theta_1 = 1.2$ (as before)
- CM3: $\theta_1 = -1.2$ (new)
- CM4: $\theta_1 = -0.3$ (new)
 
START rerun simulation: 05/06/2024

FINISH rerun simulation: 07/06/2024

### Real data - EGG

I noticed some error in my SNP annotation: when I checked the SNPs in script *01_Prep_02*, I sort by p-value to get the lead SNP per locus. But I forgot to reorder by chromosome and position. Hence, all the rsIDs were wrong in the GX and GY tables. I decided to rerun everything. 

In addition, I added the 2-sample MVMR analyses using EGG and UKBB as 2nd sample for the outcome associations: 

- **EGG consortium**: birth weight Z-scores, problematic as these are my candidate SNPs
- **UKB raw**: downloaded from Neale lab, self-reported birth weight in kg
- **UKB irnt**: downloaded from Neale lab, self-reported birth weight, inverse-rank normal transformed (quantile-like)
- **Sakaue et al.**: Cesarian section, GCST90018810 in the GWAS Catalog, Meta-analysis of UKB, FinnGen, and BioBank Japan
- **UKB elective**: downloaded from Neale lab, delivery mode of own children being elective CS (against all other modes)
- **UKB emergency**: downloaded from Neale lab, delivery mode of own children being emergency CS (against all other modes)

**Conclusion**:

- affirmative for BW: regardless of what data set was used, we can replicate the POPS findings (EFW mean and slope have positive causal effect on BW)
- negative for CS


## 2024/05/31

Okay, what did I do in the last week? 

### Simulation

I run the script to evaluate the main simulation and added some extra routines to create some nicer plots for my talk at the MR Conference in Bristol in June. 

I also added and run all sensitivity simulation checks. This includes one for weak instruments, using MV GMM instead of MV IVW. 

Because I feel quite good now about the simulation, I removed all the old files to the rds (GAMLSS out files are quite large, and are taking too much space away). This includes: 

- simulation
- simulation2
- simulation_v2
- simulation_v3
- simulation_v4

I only keep main, sensitivity and test (not tracked). 

### Real data

Everything is updated: I repeated all analyses with 100 independent SNPs from the EGG consortium data. I think about adding here the 2-sample approach. Maybe next week. 

## 2024/05/24

Yes, it has been again too long since I tracked anything here. 

### Simulation

I decided to restructure my simulation. In the **main simulation**, I want to focus on GAMLSS as longitudinal GWAS model of the exposure, real age as time parameter, quadratic growth, negative correlation between the SNP effects of the mean and the slope, with one shared SNP set for mean and slope, and a distinct SNP set for the variability effect. To have good power, I simulate 10,000 samples and 15 time points. 

I have decided on 12 scenarios: 

| Scenario | AS used to simulate $X$   | $\theta_1$ | $\theta_2$ | $\theta_3$ |
| :------- | :------------------------ | :--------- | :--------- | :--------- |
| 1        | $AS_1$ and $AS_2$         | 0.3        | 0.3        | 0          |
| 2        | $AS_1$ and $AS_2$         | 1.2        | 0.3        | 0          |
| 3        | $AS_1$ and $AS_2$         | 1.2        | 0.3        | 1          |
| 4        | $AS_1$ and $AS_2$         | 1.2        | 0          | 1          |
| 5        | $AS_1$, $AS_2$ and $AS_3$ | 0.3        | 0.3        | 0          |
| 6        | $AS_1$, $AS_2$ and $AS_3$ | 1.2        | 0.3        | 0          |
| 7        | $AS_1$, $AS_2$ and $AS_3$ | 1.2        | 0.3        | 1          |
| 8        | $AS_1$, $AS_2$ and $AS_3$ | 1.2        | 0          | 1          |
| 9        | $AS_1$ and $AS_3$         | 0.3        | 0.3        | 0          |
| 10       | $AS_1$ and $AS_3$         | 1.2        | 0.3        | 0          |
| 11       | $AS_1$ and $AS_3$         | 1.2        | 0.3        | 1          |
| 12       | $AS_1$ and $AS_3$         | 1.2        | 0          | 1          |

Later, I will run some sensitivity simulations:

1) positive correlation between the mean and slope effect
2) no correlation between the mean and slope effect
3) using GAMLSS but estimating only mean and slope effect
4) using GAMLSS but estimating only mean and variability effect
5) assuming one shared SNP set for all mean, slope and variability
6) assuming three distinct SNP sets for mean, slope and variability
7) using visit as time parameter
8) reducing sample size and time points to POPS size (3,000 samples, 5 time points)
9) using binary outcomes instead of continuous ones
10) using mvmr_gmm to correct for weak instrument bias

### Real data

I also started on rerunning the real data analysis, using the genome-wide significant results from the EGG consortium. This results in 100 candidate SNPs, which might make writing-up a bit more simpler. 

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