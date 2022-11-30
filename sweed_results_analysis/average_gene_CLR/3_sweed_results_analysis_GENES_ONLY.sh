### To run this you need: your species genes .gff file, chrom.txt file, orthogroups map.txt file and all SweeD_Report files for each species
### Of course, you also need the Rscripts called, which you should have downloaded together with this script




### This is to take care of badly formatted gff with spaces/tabs in the ID field, and adding 500bp flanks each side of each gene

awk '{OFS="\t"}{print $1,$2,$3,$4-500,$5+500,$6,$7,$8,"id=gene;"}' genes.gff > noid_genes.gff


### Taking care of genes whose start goes below 0 after adding 500bp flanks

awk '{OFS="\t"}{if ($4 > 0)     print $1,$2,$3,$4,$5,$6,$7,$8,$9;else print $1,$4+500,$3,1,$5,$6,$7,$8,$9;}' noid_genes.gff > tmp && mv tmp noid_genes.gff


### First I format SweeD outputs and add the chrom/scaffold name to all files, while removing first 3 lines (empty line, garbage and header)
### CHANGE PREFIX OF SWEED REPORT FILE FOR EACH SPECIES BEFORE RUNNING THIS SCRIPT

while read p; do

tail -n +4 SweeD_Report.mitchell_bstricta.vcf.gz_$p| awk -v var="$p" 'BEGIN {OFS="\t"} {print var, $1, $1, "CLRSCAN", $2}' > formatted\_$p

done < chrom.txt




### Concatenate formatted outputs into 1 and remove all the formatted chrom/scaffolds
cat formatted_* > all.bed

rm formatted*

################################### ADDITIONAL LINE TO CALCULATE EMP P VALUES BASED ON SCANS WITHIN GENES ONLY
/data/programs/bedtools2/bin/bedtools intersect -a all.bed -b noid_genes.gff -wa | uniq > tmp
mv tmp all.bed
######################################################################################


### Calculating empirical p values for all clr scans, including all genome-wide scans 

Rscript emp.R




### Here I make a list of genes intersecting the clr scans. Each gene gets reported as many times as the number of clr scan it includes

/data/programs/bedtools2/bin/bedtools intersect -a noid_genes.gff -b all_emp.bed -wo | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$15}' > genes_emp.gff




### Here first I add for each gene, the count of CLR scans - then I format the final file to only have gene name (given by chrom:start-end), empirical pvalue and number of CLR scans within the gene
### I also remove the 500bp flanks from each side of genes, as I need unflanked genes to match genes with orthogroups using the map I created from Jim file

/data/programs/bedtools2/bin/bedtools coverage -a genes_emp.gff -b all_emp.bed -counts | awk '{OFS="\t"}{if ($4 > 1) print $1":"$4+500"-"$5-500,$9,$10; else print $1":"$2"-"$5-500,$9,$10;}' > final_genes_analysis.txt




### I use the script below to take the average emp pvalue for each gene

Rscript average.R




### Now it is time to add the orthogroup for each gene, retrieved from map.txt

awk '{OFS="\t"}NR==FNR { id[$1]=$0; next } ($1 in id){ print id[$1], $2}' final_genes_average.txt map.txt > final_genes_average_ortho.txt




### Count genes in each orthogroup in a very patchy way

Rscript ortho_count.R

awk '{print $5}' final_genes_average_ortho.txt > ortho_list.txt

./create_file.sh > count.txt

awk '{print $2}' count.txt > final_count.txt

paste final_genes_average_ortho.txt final_count.txt > tmp_file.txt




### Add header to avoid loosing track of what each column is

sed -i '1s/^/gene\tmin_emp_p\tscan_n\tmean_emp_p\torthogroup\tortho_size\n/' tmp_file.txt




### Apply dunn sidak for each orthogroup with R: taking the gene in each orthogroup with minimum dunn_sidak adjusted emp_p value (previously calculated for each gene, corrected for clr scan number in a gene)
### Then I further correct this with dunn sidak again, but for the number of genes in the orthogroup where the gene belong
### Final file produced: final_orthogroup.txt

Rscript ortho_dunn_sidak.R




### clean the huge amount of garbage generated


rm count.txt
rm final_genes_analysis.txt
rm ortho_list.txt
rm final_count.txt
rm genes_emp.gff
rm all*
rm final_genes_average.txt
rm noid_genes.gff
rm ortho_gene_count.txt
rm tmp_file.txt
rm final_genes_average_ortho.txt
