#!/bin/bash
#SBATCH --job-name=wza_run10
#SBATCH --time=2:00:00  # Time limit set to 4 days
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30gb
#SBATCH --output=wza_run10_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

module load python/3.11_conda
conda activate /home/tbellagio/miniforge3/envs/pipeline_snakemake
export LD_LIBRARY_PATH="/home/tbellagio/miniforge3/envs/run_baypass/lib:$LD_LIBRARY_PATH"
cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/wza
python /carnegie/nobackup/scratch/tbellagio/gea_grene-net/wza/general_WZA_script_mod.py --correlations kendall_10_w_id_n_blocks.csv --summary_stat K_tau_p --window "block" --output wza_run10_results.csv --sep ","

