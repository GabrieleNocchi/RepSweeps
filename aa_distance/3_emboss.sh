#!/bin/bash
#SBATCH --time=0-01:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=500M
#SBATCH --account=def-yeaman
#SBATCH --array=1-44

FASTA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list.txt)
OUTPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list2.txt)


module load emboss
distmat -sequence $FASTA -protmethod 1 -outfile $OUTPUT\_distances.txt