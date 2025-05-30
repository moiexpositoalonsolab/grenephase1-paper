#!/bin/bash
#SBATCH --job-name=chain_2
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=25
#SBATCH --cpus-per-task=8
#SBATCH --output=chain_2_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

export PATH="${PATH}:/home/tbellagio/bin"

cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass/neutral_runs

baypass -gfile /carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass/allele_counts_3rdgen_nheader.txt -poolsizefile /carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass/pool_sizes_3rdgen_nheader.txt -nthreads 8 -seed 54894774 -print_omega_samples -outprefix chain_3


