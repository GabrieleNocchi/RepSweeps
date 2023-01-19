for FILE in *.vcf.gz
do
	#extract number of samples to calculate maf threshold to use to exclude singletons snps
	samples=$(/data/programs/bcftools-1.9/bcftools query -l $FILE| wc -l)
	af=$(bc <<< "scale=4; 2/($samples*2)")


	#filter
	/data/programs/vcftools_0.1.13/bin/vcftools --gzvcf $FILE --max-missing 0.7 --maf $af --minQ 30 --minGQ 20 --minDP 5 --max-alleles 2 --recode --recode-INFO-all


	rm out.log
	bgzip out.recode.vcf
	mv out.recode.vcf.gz $FILE
done
