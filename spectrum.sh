for FILE in *.vcf.gz
do

/data/programs/vcftools_0.1.13/bin/vcftools --counts2  --gzvcf $FILE --stdout | awk 'NR<=1 {next} {print $2"\t"$6"\t"$4"\t1"}'| sort -nk1 |awk 'BEGIN{print "position\tx\tn\tfolded"};{print $1"\t"$2"\t"$3"\t"$4}' > $FILE\_tmp.txt

/data/programs/SweeD_v3.2.1_Linux/SweeD-P -input $FILE\_tmp.txt -osfs $FILE\_spectrum.txt -strictPolymorphic -threads 8

done


wait
echo "Analysis completed"
