#!/bin/bash
#SBATCH --job-name=changep_time_28_16
#SBATCH --time=24:00:00  # 
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30gb
#SBATCH --output=changep_time_28_16_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

module load python/3.11_conda
conda activate /home/tbellagio/miniforge3/envs/pipeline_snakemake
cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/snp_origin
python calc_delta_time.py 28 16

