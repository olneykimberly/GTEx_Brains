#!/bin/bash

#SBATCH -p phi                                                                 
#SBATCH -n 28                          # Max is 256 on phi node                
#SBATCH -t 48:00:00                                                           
#SBATCH -q wildfire    
#SBATCH --job-name=PCA_T_vs_T          # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH --mail-user=mjpete11@asu.edu   # send-to address

module load r/3.5.2

Rscript PCA_Tissue.R 
