#!/bin/bash

#SBATCH --job-name=launch_simulation
#SBATCH --A wheat_striperust
#SBATCH --output=simulation_output.log
#SBATCH --error=simulation_error.log
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=compute