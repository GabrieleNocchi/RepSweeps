#!/bin/bash
#SBATCH --time=0-24:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --account=def-yeaman

module load r/4.2.1
module load bedtools
./3_sweed_results_analysis_genes_only_SweeD.sh
