#!/bin/bash
#
#SBATCH --job-name=SMARTSLAB-KT-55555deg
#SBATCH --output=5deg.log
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=2048
#SBATCH --mail-user=tienyiah@uci.edu
#
#SBATCH --array=1-73

srun julia fit_KT_fixedTd_ForwardDiff_single_longitude.jl $SLURM_ARRAY_TASK_ID
