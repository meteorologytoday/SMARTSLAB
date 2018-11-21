#!/bin/bash
#
#SBATCH --job-name=SMARTSLAB-KT-2deg
#SBATCH --output=2deg.log
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=2048
#SBATCH --mail-user=tienyiah@uci.edu

#
##SBATCH --array=1-1
#SBATCH --array=1-5

srun julia stanfit_KT_single_longitude.jl $SLURM_ARRAY_TASK_ID
