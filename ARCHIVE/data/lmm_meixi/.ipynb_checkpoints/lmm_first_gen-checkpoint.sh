#!/bin/bash
#SBATCH --job-name=lmm_first_gen
#SBATCH --time=3-00:00:00
#SBATCH --ntasks=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=21  
#SBATCH --output=result-%j.out       
#SBATCH --error=error-%j.log         
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

echo start
source /home/tbellagio/miniforge3/etc/profile.d/conda.sh

export PATH="${PATH}:/home/username/bin"

cd '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/lmm_meixi/'

conda activate /home/tbellagio/miniforge3/envs/r-environment
Rscript step1_run_lmm.R
