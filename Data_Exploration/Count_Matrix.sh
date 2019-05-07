#!/bin/bash
#SBATCH --job-name=Count_Matrix # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH --mail-user=mjpete11@asu.edu   # send-to address
#SBATCH -p phi
#SBATCH -n 150                         # Max is 256 on phi node
#SBATCH -t 48:00:00
#SBATCH --qos=normal

module load r/3.5.2

Rscript Count_Matrix.R --keep-target-files --rerun-incomplete --cluster "sbatch -p phi -n 150 -t 48:00:00"
