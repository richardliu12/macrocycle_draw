#!/bin/bash

#SBATCH -J csearch

#SBATCH -n 8           # Number of cores 

#SBATCH -N 1           # 4 cores per node

#SBATCH -p jacobsen 

#SBATCH --mem 4000      

#SBATCH -t 2-00:00      # Two day max job

#SBATCH -e error.err    # Standard error file

/n/sw/schrodinger/macromodel -WAIT csearch"${SLURM_ARRAY_TASK_ID}".com
