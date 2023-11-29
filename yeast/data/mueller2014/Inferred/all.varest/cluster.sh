#!/bin/bash

#SBATCH --job-name=ann-emp
#SBATCH --partition compute
#SBATCH --time=4-0
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=anneal.log
#SBATCH --error=anneal.log

python3 code.py 3 2e2 5e-10 2000000 5000 empirical
