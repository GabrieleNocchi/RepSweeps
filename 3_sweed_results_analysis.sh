### need genes .gff file, chrom.txt file, orthogroups map.txt file and all SweeD_Report files for each species

# This is to take care of badly formatted gff with spaces/tabs in the ID field, and adding 500bp flanks each side of each gene

awk '{OFS="\t"}{print $1,$2,$3,$4-500,$5+500,$6,$7,$8,"id=gene;"}' genes.gff > noid_genes.gff

# first I format SweeD outputs and add chrom/scaffold name to all files, while removing first 3 lines (empty line, garbage, header)
while read p; do

tail -n +4 SweeD_Report.lowry_phallii.vcf.gz_$p| awk -v var="$p" 'BEGIN {OFS="\t"} {print var, $1, $1, "CLRSCAN", $2}' > formatted\_$p

done < chrom.txt

#concatenate formatted outputs into 1 and remove all the formatted chrom/scaffolds
cat formatted_* > all.bed

rm formatted*


#calculating empirical p values for all clr scans
Rscript emp.R


# here we make a list of genes intersecting the clr scans. Each gene gets reported as many times as the number of clr scan it includes

/data/programs/bedtools2/bin/bedtools intersect -a noid_genes.gff -b all_emp.bed -wo | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$15}' > genes_emp.gff


# here first I add for each gene, the count of CLR scans - then I format the final file to only have gene name (given by chrom:start-end), empirical pvalue and number of CLR scans within the gene
# also removing the 500bp flanks from each side of genes, as I need unflanked genes to match genes with orthogroups using the map we have

/data/programs/bedtools2/bin/bedtools coverage -a genes_emp.gff -b all_emp.bed -counts | awk '{OFS="\t"}{print $1":"$4+500"-"$5-500,$9,$10}' > final_genes_analysis.txt


# I use the script below to take the minimum emp pvalue for each gene and apply a dunn sidak corrrection to it, based on the nummber of CLR scan within each gene
# final file produced = final_genes_dunak.txt
Rscript dunn_sidak.R


#now I add the orthogroup for each gene, retrieved from map.txt
awk '{OFS="\t"}NR==FNR { id[$1]=$0; next } ($1 in id){ print id[$1], $2}' final_genes_dunak.txt map.txt > final_genes_dunak_ortho.txt
# add header
sed -i '1s/^/gene\tmin_emp_p\tscan_n\tdunn_sidak\torthogroup\n/' final_genes_dunak_ortho.txt

# clean a bit
rm genes_emp.gff
rm all*
rm final_genes_analysis.txt
rm noid_genes.gff
