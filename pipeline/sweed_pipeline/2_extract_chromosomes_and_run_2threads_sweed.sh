#First I create a neutral SFS for each dataset, based on all scaffolds/chromosomes
#then for each vcf, I make separate chromosomes/scaffolds vcf;
#then on each chromosome/scaffold SweeD is run using the grid specified in the grid file previously generated and the neutral SFS based on the entire genome

regex='(.+)	(.+)'

for FILE in *.vcf.gz
do

/data/programs/vcftools_0.1.13/bin/vcftools --counts2  --gzvcf $FILE --stdout | awk 'NR<=1 {next} {print $2"\t"$6"\t"$4"\t1"}'| sort -nk1 |awk 'BEGIN{print "position\tx\tn\tfolded"};{print $1"\t"$2"\t"$3"\t"$4}' > $FILE\_tmp.txt

/data/programs/SweeD_v3.2.1_Linux/SweeD-P -input $FILE\_tmp.txt -osfs $FILE\_spectrum.txt -strictPolymorphic -threads 8


 while read p; do
                

 	if [[ $p =~ $regex ]];  then    /data/programs/bcftools-1.9/bcftools view $FILE --regions ${BASH_REMATCH[1]} -Ov -o $FILE\_${BASH_REMATCH[1]}\.vcf; /data/programs/SweeD_v3.2.1_Linux/SweeD-P -input $FILE\_${BASH_REMATCH[1]}\.vcf -folded -grid ${BASH_REMATCH[2]} -name $FILE\_${BASH_REMATCH[1]} -strictPolymorphic -isfs $FILE\_spectrum.txt -threads 2; fi
 done < $FILE\_genome_grid_1k.txt


done
wait
echo "Analysis completed"
