for FILE in *.vcf.gz
do

/data/programs/vcftools_0.1.13/bin/vcftools --counts2  --gzvcf $FILE --stdout | awk 'NR<=1 {next} {print $2"\t"$6"\t"$4"\t1"}'| sort -nk1 |awk '{print $1"\t"$2"\t"$3"\t"$4}'> $FILE\_tmp.txt

shuf -n 100000 *tmp* > tmp

sort -nk1 tmp > tmp2 
sed -i '1s/^/position\tx\tn\tfolded\n/' tmp2 
mv tmp2 $FILE\_tmp.txt
rm tmp



/data/programs/SweeD_v3.2.1_Linux/SweeD-P -input $FILE\_tmp.txt -osfs spectrum.txt -strictPolymorphic -threads 2

done


wait
echo "Analysis completed"
