#!/bin/bash
#SBATCH --partition=stat
#SBATCH --time=2:00:00
#SBATCH --nodes=1 --ntasks-per-node=25 --mem=5gb
#SBATCH --output=/home/zhang.7824/Meijia-two-sample-test-Revision-1/pbs_logs/new-subsample-alpha175-spar050-dgn1-nbm1-log.txt
#SBATCH --job-name=subsample-aa175-050-dgn1-nbm1
#SBATCH --mail-user=ronaldaylmerfisher@gmail.com

#SBATCH --mail-type=ALL

module load matlab

cd ~/Meijia-two-sample-test-Revision-1/
matlab -r "alpha=1.75;GraphonName1='NewDegenGraphon1';GraphonName2='NewBlockModel1';A_sparse_power=-0.50;B_sparse_power=-0.50;new_degeneracy_no_alpha_benchmark_subsample"
