#!/bin/bash
#SBATCH --job-name=wza_run5
#SBATCH --time=2:00:00  # 
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30gb
#SBATCH --output=wza_run5_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

module load python/3.11_conda
conda activate /home/tbellagio/miniforge3/envs/pipeline_snakemake
export LD_LIBRARY_PATH="/home/tbellagio/miniforge3/envs/run_baypass/lib:$LD_LIBRARY_PATH"
cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/linages_wza_picmin
python /carnegie/nobackup/scratch/tbellagio/gea_grene-net/wza/general_WZA_script_mod.py --correlations kendall_5_w_id_n_blocks.csv --summary_stat K_tau_p --window "block" --output wza_run_5_results.csv --sep ","

