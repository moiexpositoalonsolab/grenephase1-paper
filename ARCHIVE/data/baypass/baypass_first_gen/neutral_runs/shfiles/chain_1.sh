#!/bin/bash
#SBATCH --job-name=chain_1
#SBATCH --time=2-00:00:00  # Time limit set to 4 days
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30gb
#SBATCH --cpus-per-task=4
#SBATCH --output=chain_1_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

module load python/3.11_conda
conda activate /home/tbellagio/miniforge3/envs/run_baypass
export LD_LIBRARY_PATH="/home/tbellagio/miniforge3/envs/run_baypass/lib:$LD_LIBRARY_PATH"
cd /home/tbellagio/scratch/gea_grene-net/baypass_terminal/neutral_runs
/home/tbellagio/bin/baypass -gfile /carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/neutral_runs/allele_counts_1stgen_nheader.txt -seed 70069661 -nthreads 4 -print_omega_samples -outprefix chain_1


