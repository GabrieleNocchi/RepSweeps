#!/bin/bash
#SBATCH --time=00-20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=500M
#SBATCH --account=def-yeaman
#SBATCH --array=1-999



GENEFILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" genes_files_list.txt)

module load bcftools





#change file name accordingly
FILE="kubota_ahalleri.vcf.gz"


regex='(.+)	(.+)'

while read p; do

if [[ $p =~ $regex ]];  then 	bcftools view $FILE --regions ${BASH_REMATCH[1]} -Ov -o $FILE\_${BASH_REMATCH[2]}\.vcf; /home/gabnoc/projects/def-yeaman/gabnoc/OmegaPlus/OmegaPlus-M -input $FILE\_${BASH_REMATCH[2]}\.vcf -minwin 200 -maxwin 100000 -grid 5 -name $FILE\_${BASH_REMATCH[2]}  -seed 12345 -threads 2; fi
	rm $FILE\_${BASH_REMATCH[2]}\.vcf
done < $GENEFILE


