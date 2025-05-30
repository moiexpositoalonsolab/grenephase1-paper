#!/bin/bash
#SBATCH --job-name=picmin
#SBATCH --time=1-00:00:00 
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30gb
#SBATCH --output=picmin_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

source /home/tbellagio/miniforge3/etc/profile.d/conda.sh

export PATH="${PATH}:/home/username/bin"

cd '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/picmin'

conda activate /home/tbellagio/miniforge3/envs/r-environment
Rscript run_pic_min.r


