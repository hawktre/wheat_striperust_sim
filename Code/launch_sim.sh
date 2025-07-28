#!/bin/bash -l

#SBATCH -J wheat_striperust_sim
#SBATCH -o wheat_striperust.output
#SBATCH -e error.output

# Default in slurm
#SBATCH --mail-user hawkintr@oregonstate.edu
#SBATCH --mail-type=ALL

# Request 4 hours run time
#SBATCH -t 24:0:0

# Specify the partition to run on
#SBATCH --partition=dri.q

# Set Memory for job
#SBATCH --mem=16G
#SBATCH --cpus-per-task=32
#SBATCH --ntasks=1

# Change version of R

echo "start R job"

module load R/4.4.1

Rscript --vanilla Code/04b_RunSim.R

echo "R Finished"