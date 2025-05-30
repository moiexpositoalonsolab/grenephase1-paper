#!/bin/bash
#SBATCH --job-name=pic_min_run100000
#SBATCH --time=2-00:00:00  # Time limit set to 4 days
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30gb
#SBATCH --output=pic_min_run100000_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

source /home/tbellagio/miniforge3/etc/profile.d/conda.sh
export PATH="${PATH}:/home/username/bin"
cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/linages_wza_picmin
conda activate /home/tbellagio/miniforge3/envs/r-environment
Rscript run_pic_min_deltap.r 100000

