#!/bin/bash
#SBATCH --job-name=split_34
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=120gb
#SBATCH --cpus-per-task=2
#SBATCH --output=split_34_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

source /home/tbellagio/miniforge3/etc/profile.d/conda.sh

export PATH="${PATH}:/home/username/bin"

cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/leave_1_out

conda activate /home/tbellagio/miniforge3/envs/r-environment
Rscript run_lfmm_full_genome_ridge_deltapmean.r 34 "['28' '32' '52' '49' '33' '60' '46' '55' '26' '10' '23' '42' '43' '37'
 '57' '6' '48' '27' '11' '45' '53' '2' '4' '12' '24']"
