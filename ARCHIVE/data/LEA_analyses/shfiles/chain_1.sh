#!/bin/bash
#SBATCH --job-name=chain_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=100gb
#SBATCH --cpus-per-task=10
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/gea_grene-net/LEA_analyses/shfiles/chain_1_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

module load python/3.11_conda
conda activate /home/tbellagio/miniforge3/envs/snakemake_pipeline
cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/LEA_analyses
Rscript run_lfmm.R allele_freq_part-00.csv /carnegie/nobackup/scratch/tbellagio/gea_grene-net/LEA_analyses/full/result_allele_freq_part-00.out
