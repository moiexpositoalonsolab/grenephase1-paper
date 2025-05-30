#!/bin/bash
#SBATCH --job-name=run_partition_binomial_reg154_30
#SBATCH --time=1:00:00  # Time limit set to 4 hours
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30gb
#SBATCH --output=run_partition_binomial_reg_154_30_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

module load python/3.11_conda
conda activate /home/tbellagio/miniforge3/envs/pipeline_snakemake
export LD_LIBRARY_PATH="/home/tbellagio/miniforge3/envs/run_baypass/lib:$LD_LIBRARY_PATH"
cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/binomial_regression_firstgen_go
python run_partition_binomial_reg_first_gen.py 154 30

