gzip -d *.vcf.gz

#for each vcf I create a genome grid file with 2 columns
#the first column has chromosomes/scaffolds names
#the second column has the grid numbers for SweeD 
#I want 1 measurement every 1kb
# First scan is at 1st SNP, last scan is at last SNP
# So I take the length between last and first SNP in each scaffold/chromosomes, and divide that by 1000


for FILE in *.vcf
do
    ./grid.pl $FILE | awk '{OFS="\t"}{print $1,int(($3-$2)/1000)}' > $FILE\.gz_genome_grid_1k.txt
done



#here I take care of exceptions, scaffolds or contigs shorter than 1k, I want to scan at least 2 positions (first and last SNP)
#the minimum grid number possible in SweeD is 2, the first scan is taken at the first SNP and the last scan is the last SNP, hence you need an extra grid (+1) to cover evenly every 1kb


for FILE2 in *_genome_grid_1k.txt
do
	awk '{OFS="\t"}{if ($2 >= 1)	print $1,$2+1;else print $1,2;}' $FILE2 >tmp && mv tmp $FILE2
done



for FILE3 in *.vcf
do
bgzip $FILE3
tabix -p vcf $FILE3\.gz
done
