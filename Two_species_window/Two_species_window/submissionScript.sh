#!/bin/bash

#SBATCH --clusters=all
#SBATCH --partition=short
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1 
#SBATCH --mem-per-cpu=8G 
#SBATCH --job-name=varying_window_RUN
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francesca.lovell-read@merton.ox.ac.uk
#SBATCH --output=out_ex.out
#SBATCH --array=1-289

module load MATLAB/R2022a #%load MATLAB (can specify version I think if need be)
matlab -nodisplay -nosplash -r "varying_window_RUN($SLURM_ARRAY_TASK_ID); quit;"