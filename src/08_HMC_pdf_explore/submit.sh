#!/bin/bash
#
#SBATCH --job-name=SMARTSLAB-KT-2deg
#SBATCH --output=2deg.log
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=2048
#SBATCH --mail-user=tienyiah@uci.edu

#
#SBATCH --array=1-73

srun julia stanfit_KT_pdf_explore.jl $SLURM_ARRAY_TASK_ID 36
