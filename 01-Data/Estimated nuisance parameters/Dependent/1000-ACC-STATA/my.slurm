#!/bin/bash

#SBATCH --mail-type=ALL  # Alerts sent when job begins, ends, or aborts
#SBATCH --nodes=1   # Number of nodes required
#SBATCH --ntasks=1   # Number of nodes required
#SBATCH --mem=2G  # Total memory (RAM) required, per node
#SBATCH --time=0-00:5:00  # Wall Clock time (dd-hh:mm:ss) [max of 14 days]
#SBATCH --array=0-1000
#SBATCH --job-name=1000-ACC # job name

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load stata
stata main2.do $SLURM_ARRAY_TASK_ID
