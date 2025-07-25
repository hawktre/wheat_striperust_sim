#!/bin/bash -l

#SBATCH -J my_job
#SBATCH -o my_job.output
#SBATCH -e error.output

# Default in slurm
#SBATCH -D ./
#SBATCH --mail-user hawkintr@oregonstate.edu
#SBATCH --mail-type=ALL

# Request 4 hours run time
#SBATCH -t 24:0:0

# Specify the project for job

#SBATCH -A wheat_striperust_sim

# Set Memory for job
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1

# Change version of R

module load R/4.4.1

echo "start R job"

module load R

Rscript --vanilla Code/04b_RunSim.R

echo “R Finished"