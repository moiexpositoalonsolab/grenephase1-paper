#!/bin/bash
#SBATCH --job-name=chain_5
#SBATCH --time=2-00:00:00  # Time limit set to 4 days
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30gb
#SBATCH --cpus-per-task=4
#SBATCH --output=chain_5_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

module load python/3.11_conda
conda activate /home/tbellagio/miniforge3/envs/run_baypass
export LD_LIBRARY_PATH="/home/tbellagio/miniforge3/envs/run_baypass/lib:$LD_LIBRARY_PATH"
cd /home/tbellagio/scratch/gea_grene-net/baypass_terminal/neutral_runs_terminal_gen
/home/tbellagio/bin/baypass -gfile /carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass_terminal/neutral_runs_terminal_gen/allele_counts_lastgen_nheader.txt -seed 42949660 -nthreads 4 -print_omega_samples -outprefix chain_5


