#!/bin/bash

## Example SLURM script for BSU cclake jobs

## Section 1: SLURM Commands

## All SLURM commands must be placed at the start of the file
## Full documentation can be found here: https://slurm.schedmd.com/sbatch.html

## Enter a short name for the job, to be shown in SLURM output
#SBATCH -J LongMR_Sim_eval

## Enter the wall-clock time limit for for your jobs.
## If jobs reach this limit they are automatically killed.
## Maximum value 36:00:00.
#SBATCH --time=00:10:00

## For single-core jobs, this number should be '1'. 
## If your job has built-in parallelism, eg using OpenMP or 
## R's foreach() and doParallel(), increase this number as desired.
## The maximum value is 56.
#SBATCH --cpus-per-task=1

## Each task is allocated 3.4G (cclake) or 6.8G (cclake-himem). 
## If this is insufficient, uncomment and edit this line.
## Maximum value 191360M or 384960M.
## #SBATCH --mem=24G

## The system can send emails when your job starts and stops.
## Values include BEGIN, END, ALL, and TIME_LIMIT_80 and TIME_LIMIT_90 
## (reaching 80% or 90% of time limit.) Specify ARRAY_TASKS to receive
## a separate mail for each task. Multiple values can be given, separated by a comma.
#SBATCH --mail-type=FAIL

## The project account name.
## Use mrc-bsu-sl2-cpu for cclake and mrc-bsu-sl2-gpu for pascal
#SBATCH -A mrc-bsu-sl2-cpu

## The partition. Use cclake for normal jobs, or cclake-himem if needed.
#SBATCH -p cclake

## GPU jobs only:
## Uncomment and specify the number of GPUs required per node, maximum 4.
## Note that there is a maximum of 3 cores per GPU.
## #SBATCH --gres=gpu:4

## Array jobs:
## #SBATCH --array=1-16

##  - - - - - - - - - - - - - -

## Section 2: Modules

# All scripts should include the first three lines.

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

# Load the latest R version.
module load R/4.2.2 gcc/9

start_time=`date +%H`

R CMD BATCH --vanilla ../scripts/Eval_Simulation.R ../scripts/Eval_Simulation.R.out
cp Rplots.pdf Eval_Simulation_RPlots.pdf
rm Rplots.pdf

R CMD BATCH --vanilla ../scripts/MT_MainTable1_SimResults.R ../scripts/MT_MainTable1_SimResults.R.out

R CMD BATCH --vanilla ../scripts/MF_MainFigure1_Power.R ../scripts/MF_MainFigure1_Power.R.out
cp Rplots.pdf MF_MainFigure1_Power.pdf
rm Rplots.pdf

R CMD BATCH --vanilla ../scripts/MF_MainFigure2_Bias.R ../scripts/MF_MainFigure2_Bias.R.out
cp Rplots.pdf MF_MainFigure2_Bias.pdf
rm Rplots.pdf

end_time=`date +%H`
echo execution time total was `expr $end_time - $start_time` hours.
