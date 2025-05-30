#!/bin/bash
#SBATCH --job-name=gemma_bio1
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=120gb
#SBATCH --cpus-per-task=2
#SBATCH --output=gemma_bio1_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

source /home/tbellagio/miniforge3/etc/profile.d/conda.sh

export PATH="${PATH}:/home/username/bin"

cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/allele_assoc_runs/only_grenenet/lmm_gemma/bio1/

conda activate /home/tbellagio/miniforge3/envs/gwas

gemma -bfile greneNet_final_v1.1.recode -maf 0.05 -lmm -k /carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/output/greneNet_final_v1.1.recode_maf05.cXX.txt -c /carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/allele_assoc_runs/only_grenenet/pc_300.txt -o bio1

