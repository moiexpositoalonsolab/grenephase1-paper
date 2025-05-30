#!/bin/bash
#SBATCH --job-name=split_1
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=120gb
#SBATCH --cpus-per-task=2
#SBATCH --output=split_1_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

source /home/tbellagio/miniforge3/etc/profile.d/conda.sh

export PATH="${PATH}:/home/username/bin"

cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/leave_1_out

conda activate /home/tbellagio/miniforge3/envs/r-environment
Rscript run_lfmm_full_genome_ridge_deltapmean.r 1 "['1' '54' '57' '27' '5' '9' '45' '23' '4' '53' '55' '43' '28' '52' '26'
 '6' '10' '32' '49' '24' '42' '46' '25' '60' '13']"
