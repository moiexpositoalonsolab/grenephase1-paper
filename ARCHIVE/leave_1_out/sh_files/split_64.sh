#!/bin/bash
#SBATCH --job-name=split_64
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=120gb
#SBATCH --cpus-per-task=2
#SBATCH --output=split_64_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

source /home/tbellagio/miniforge3/etc/profile.d/conda.sh

export PATH="${PATH}:/home/username/bin"

cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/leave_1_out

conda activate /home/tbellagio/miniforge3/envs/r-environment
Rscript run_lfmm_full_genome_ridge_deltapmean.r 64 "['54' '53' '2' '1' '49' '60' '28' '45' '5' '10' '13' '24' '11' '42' '43'
 '23' '4' '48' '12' '27' '9' '6' '46' '57' '25']"
