
regex='(.+)	(.+)'

for FILE in *.vcf.gz
do


while read p; do
                

 	if [[ $p =~ $regex ]];  then    /data/programs/bcftools-1.9/bcftools view $FILE --regions ${BASH_REMATCH[1]} -Ov -o $FILE\_${BASH_REMATCH[1]}\.vcf; /data/programs/omegaplus-master/OmegaPlus-M -input $FILE\_${BASH_REMATCH[1]}\.vcf -minwin 200 -maxwin 10000 -grid ${BASH_REMATCH[2]} -name $FILE\_${BASH_REMATCH[1]}  -seed 12345 -threads 2; fi
done < $FILE\_genome_grid_1k.txt


done
wait
echo "Analysis completed"
