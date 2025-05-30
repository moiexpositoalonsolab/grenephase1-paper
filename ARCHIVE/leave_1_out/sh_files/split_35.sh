#!/bin/bash
#SBATCH --job-name=split_35
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=120gb
#SBATCH --cpus-per-task=2
#SBATCH --output=split_35_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

source /home/tbellagio/miniforge3/etc/profile.d/conda.sh

export PATH="${PATH}:/home/username/bin"

cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/leave_1_out

conda activate /home/tbellagio/miniforge3/envs/r-environment
Rscript run_lfmm_full_genome_ridge_deltapmean.r 35 "['2' '55' '45' '57' '24' '25' '13' '27' '37' '23' '32' '54' '46' '33' '48'
 '26' '53' '43' '9' '11' '60' '12' '6' '1' '42']"
