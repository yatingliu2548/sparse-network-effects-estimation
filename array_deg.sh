#!/bin/bash

#SBATCH --job-name=sim
#SBATCH --output=experiment/run_sim_%A_%a.out
#SBATCH --error=experiment/run_sim_%A_%a.err
#SBATCH --array=1-100
#SBATCH --time=24:00:00
#SBATCH --partition=cdonnat
#SBATCH --ntasks=1
#SBATCH --mem=6G
#SBATCH --account=pi-cdonnat
#SBATCH --qos=cdonnat
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yatingliu@rcc.uchicago.edur

# Print the task id.
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "My SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
# Add lines here to run your computations
job_id=$SLURM_ARRAY_JOB_ID
#module load gsl
#module load R/4.2.0
#module load matlab
#module load python
source activate ${SCRATCH}/${USER}/.local/share/r-miniconda/envs/r-reticulate
#sinteractive
module load gsl
module load matlab
module load python
module unload R
module load R/4.2.0
source activate ${SCRATCH}/${USER}/.local/share/r-miniconda/envs/r-reticulate
#source activate ~/.local/share/r-miniconda/envs/r-reticulate

MATLAB_PATH="/software/matlab-2023a-el8-x86_64/bin/matlab"
result_file="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_$1"
#result_file="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_$1_$2"
#result_file="${SLURM_ARRAY_JOB_ID}_${1234}_$1"
echo "result file is ${result_file}"
cd $SCRATCH/$USER/sparse-network-effects-estimation/
#cd topic-modeling/
working_dir="${SCRATCH}/${USER}/sparse-network-effects-estimation/"
#working_dir="topic-modeling/"
#Rscript synthetic_array.R $1234 $result_file $15 $MATLAB_PATH # 5 topic
matlab -nodisplay -r "run('eta3_array_deg.m');exit" 
# 5 topic