#!/bin/bash
#SBATCH --job-name=split_53
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=120gb
#SBATCH --cpus-per-task=2
#SBATCH --output=split_53_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

source /home/tbellagio/miniforge3/etc/profile.d/conda.sh

export PATH="${PATH}:/home/username/bin"

cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/leave_1_out

conda activate /home/tbellagio/miniforge3/envs/r-environment
Rscript run_lfmm_full_genome_ridge_deltapmean.r 53 "['26' '28' '9' '33' '43' '54' '12' '1' '42' '55' '2' '24' '11' '27' '6'
 '4' '57' '32' '52' '13' '23' '10' '5' '46' '60']"
