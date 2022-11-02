#!/bin/bash
#SBATCH --time=1-00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=500M
#SBATCH --account=def-yeaman
#SBATCH --array=1-16

#create chrom.txt and grid.txt with awk from the *_genome_grid_1k.txt file previously generated for each species
CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" chrom.txt)
GRID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" grid.txt)

module load bcftools
#change file name accordingly
FILE="wright_atuber.vcf.gz"


bcftools view $FILE --regions $CHROM -Ov -o $FILE\_$CHROM\.vcf; /home/gabnoc/projects/def-yeaman/gabnoc/SweeD/SweeD-P -input $FILE\_$CHROM\.vcf -folded -grid $GRID -name $FILE\_$CHROM -strictPolymorphic -isfs $FILE\_spectrum.txt -threads 2


