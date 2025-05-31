#!/bin/bash
#SBATCH --job-name=wzasite_54
#SBATCH --time=1:00:00  # Time limit set to 4 days
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30gb
#SBATCH --output=wza_site_54%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

module load python/3.11_conda
conda activate /home/tbellagio/miniforge3/envs/pipeline_snakemake
export LD_LIBRARY_PATH="/home/tbellagio/miniforge3/envs/run_baypass/lib:$LD_LIBRARY_PATH"
cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/gwas/fitness_gwas
python general_WZA_script_mod_polynomial_order7.py --correlations site_54/output/results_lmm.csv --summary_stat p_wald --window blocks --output site_54/wza_fitness_gwas_poly7.csv --sep ","


