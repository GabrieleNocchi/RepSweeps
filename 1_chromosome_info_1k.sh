gzip -d *.vcf.gz

#for each vcf I create a genome grid file with 2 columns
#the first column has chromosomes/scaffolds names
#the second column has the grid numbers for SweeD and Omegaplus. 
#I want 1 measurement every 10kb

for FILE in *.vcf
do
    awk -F "=|," '{OFS="\t"}($4~/length/){print $3,int($5/1000)}' $FILE | sed 's/>$//' > $FILE\.gz_genome_grid_1k.txt
done

#here I take care of exceptions, scaffolds or contigs shorter than 2k, I want to scan at least 3 positions
#the minimum grid number possible in both tools is 2, the first scan is taken at the first SNP and the last scan is the last SNP, hence I want a minimum of 3 to also have a scan in the middle
#in addition, for the chromosomes/scaffolds longer than 2k I added an additional scan (+1), as the first scan is first SNP and last scan is last SNP, so you need one more to cover
#properly every 1k bases. This is not super important, not adding it would result in one less scan at the end of the chromosomes/scaffolds.

for FILE2 in *_genome_grid_1k.txt
do
	awk '{OFS="\t"}{if ($2 > 2)	print $1,$2+1;else print $1,3;}' $FILE2 >tmp && mv tmp $FILE2
done

for FILE3 in *.vcf
do
bgzip $FILE3
tabix -p vcf $FILE3\.gz
done
