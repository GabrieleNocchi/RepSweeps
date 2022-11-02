### need genes .gff file, chrom.txt file and all SweeD_Report files for each species

# first I format SweeD outputs and add chrom/scaffold name to all files, while removing first 3 lines (empty line, garbage, header)
while read p; do

tail -n +4 SweeD_Report.wright_atuber.vcf.gz_$p| awk -v var="$p" 'BEGIN {OFS="\t"} {print var, $1, $1, "CLRSCAN", $2}' > formatted\_$p

done < chrom.txt

#concatenate formatted outputs into 1 and remove all the formatted chrom/scaffolds
cat formatted_* > all.bed

rm formatted*


#calculating empirical p values for all clr scans
Rscript emp.R


# here we make a list of genes intersecting the clr scans. Each gene gets reported as many times as the number of clr scan it includes

/data/programs/bedtools2/bin/bedtools intersect -a *.gff -b all_emp.bed -wo | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$15}' > genes_emp.gff


# here first I add for each gene, the count of CLR scans - then I format the final file to only have gene name (given by chrom:start-end), empirical pvalue and number of CLR scans within the gene

/data/programs/bedtools2/bin/bedtools coverage -a genes_emp.gff -b all_emp.bed -counts | awk '{OFS="\t"}{print $1":"$4"-"$5,$9,$10}' > final_genes_analysis.txt


# I use the script below to take the minimum emp pvalue for each gene and apply a dunn sidak corrrection to it, based on the nummber of CLR scan within each gene
# final file produced = final_genes_dunak.txt
Rscript dunn_sidak.R


# clean a bit
rm genes_emp.gff
rm all*
rm final_genes_analysis.txt
