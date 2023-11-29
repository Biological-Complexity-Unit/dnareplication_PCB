#!/bin/bash

#SBATCH --job-name=ann-emp-reps
#SBATCH --partition compute
#SBATCH --time=4-0
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks=1
#SBATCH --array=1-10
#SBATCH --cpus-per-task=1
#SBATCH --output=reps/%a/anneal.log
#SBATCH --error=reps/%a/anneal.log

python3 code.py $[30+$SLURM_ARRAY_TASK_ID-1] 2e2 5e-10 2000000 5000 empirical "reps/$SLURM_ARRAY_TASK_ID/"
