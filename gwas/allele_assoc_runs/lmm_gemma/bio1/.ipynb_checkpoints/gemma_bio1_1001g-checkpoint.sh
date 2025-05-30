#!/bin/bash
#SBATCH --job-name=gemma_bio1_1001g
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=120gb
#SBATCH --cpus-per-task=2
#SBATCH --output=gemma_bio1_1001g_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

source /home/tbellagio/miniforge3/etc/profile.d/conda.sh

export PATH="${PATH}:/home/username/bin"

cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/allele_assoc_runs/lmm_gemma/bio1/

conda activate /home/tbellagio/miniforge3/envs/gwas

gemma -bfile 1001g_grenet_climate -maf 0.05 -lmm -k /carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/allele_assoc_runs/1001g_grenet_climateLDpruned_05maf.cXX.txt -c /carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/allele_assoc_runs/pc_1000.txt -o bio1

