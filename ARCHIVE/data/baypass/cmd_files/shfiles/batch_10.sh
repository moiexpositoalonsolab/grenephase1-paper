#!/bin/bash
#SBATCH --job-name=batch_10
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=25
#SBATCH --cpus-per-task=48
#SBATCH --output=batch_10_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

export PATH="${PATH}:/home/tbellagio/bin"

cd /carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass/cmd_files/results/

baypass -gfile /carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass/individual_gfiles/partition_10.txt -efile /carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass/env_3rdgen_nheader.txt -omegafile /carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass/omegaavg_nheader.txt -poolsizefile /carnegie/nobackup/scratch/tbellagio/gea_grene-net/baypass/pool_sizes_3rdgen_nheader.txt -outprefix partition_10_chain_1 -d0yij 0.4 -seed 14077 -pilotlength 1000 -nval 50000 -npilot 40 
