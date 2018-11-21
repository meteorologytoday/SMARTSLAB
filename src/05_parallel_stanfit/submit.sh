#!/bin/bash
#
#SBATCH --job-name=SMARTSLAB-KT-5deg
#SBATCH --output=5deg.log
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=2048
#SBATCH --mail-user=tienyiah@uci.edu

#
##SBATCH --array=1-1
#SBATCH --array=2-3

srun julia stanfit_KT_single_longitude.jl $SLURM_ARRAY_TASK_ID
