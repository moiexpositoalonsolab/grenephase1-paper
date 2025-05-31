#!/bin/bash
#SBATCH --job-name=wza
#SBATCH --time=1:00:00  # Time limit set to 4 days
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30gb
#SBATCH --output=wza_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

module load python/3.11_conda
conda activate /home/tbellagio/miniforge3/envs/pipeline_snakemake
export LD_LIBRARY_PATH="/home/tbellagio/miniforge3/envs/run_baypass/lib:$LD_LIBRARY_PATH"
cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/wza_last_gen
    python general_WZA_script_mod_polynomial_order7.py         --correlations /carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/allele_assoc_runs/lmm_gemma/bio7/output/results_lmm_maf.csv         --summary_stat p_wald         --window blocks         --output gemma_1001_wza_bio7.csv --sep ","


