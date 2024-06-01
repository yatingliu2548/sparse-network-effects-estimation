#!/bin/bash
#SBATCH --partition=stat
#SBATCH --time=2:00:00
#SBATCH --nodes=1 --ntasks-per-node=25 --mem=5gb
#SBATCH --output=/home/zhang.7824/Meijia-two-sample-test-Revision-1/pbs_logs/new-resample-alpha150-spar025-dgn2-nsm2-log.txt
#SBATCH --job-name=resample-aa150-025-dgn2-nsm2
#SBATCH --mail-user=ronaldaylmerfisher@gmail.com

#SBATCH --mail-type=ALL

module load matlab

cd ~/Meijia-two-sample-test-Revision-1/
matlab -r "alpha=1.50;GraphonName1='NewDegenGraphon2';GraphonName2='NewSmoothGraphon2';A_sparse_power=-0.250;B_sparse_power=-0.250;new_degeneracy_no_alpha_benchmark_resample"