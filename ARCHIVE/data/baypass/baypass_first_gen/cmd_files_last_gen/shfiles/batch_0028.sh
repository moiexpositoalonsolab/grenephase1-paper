#!/bin/bash
#SBATCH --job-name=baypass_cmd_0028
#SBATCH --partition=expansion
#SBATCH --sockets-per-node=2
#SBATCH --cores-per-socket=32
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=25
#SBATCH --cpus-per-task=8
#SBATCH --output=batch_0028_%j.out
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

cd /central/scratch/tbellagi/gea/baypass_terminal/results
module load parallel/20220522-gcc-11.3.1-ki3yyer

cat /central/scratch/tbellagi/gea/baypass_terminal/cmd_files/run_01/catfile_0028.txt | parallel -j 8 "{} > output_catfile0028_{%}.out 2> output_catfile0028_{%}.err"

