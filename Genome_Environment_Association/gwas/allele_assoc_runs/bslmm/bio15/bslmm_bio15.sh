#!/bin/bash
#SBATCH --job-name=bslmm_bio15
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=120gb
#SBATCH --cpus-per-task=2
#SBATCH --output=bslmm_bio15_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

source /home/tbellagio/miniforge3/etc/profile.d/conda.sh

export PATH="${PATH}:/home/username/bin"

cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/idea_fromind_to_pop/allele_assoc_runs/bslmm/bio15/

conda activate /home/tbellagio/miniforge3/envs/gwas

gemma -bfile 1001g_grenet_climate -maf 0.05 -bslmm 1 -o bio15

