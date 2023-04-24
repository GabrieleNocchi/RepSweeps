### To run this you need: your species genes .gff file, chrom.txt file, orthogroups map.txt file and all SweeD_Report files for each species
### Of course, you also need the Rscripts called, which you should have downloaded together with this script



### This is to take care of badly formatted gff with spaces/tabs in the ID field, and adding 500bp flanks each side of each gene

awk '{OFS="\t"}{print $1,$2,$3,$4-500,$5+500,$6,$7,$8,"id=gene;"}' genes.gff > noid_genes.gff


### Taking care of genes whose start goes below 0 after adding 500bp flanks

awk '{OFS="\t"}{if ($4 > 0)     print $1,$2,$3,$4,$5,$6,$7,$8,$9;else print $1,$4+500,$3,1,$5,$6,$7,$8,$9;}' noid_genes.gff > tmp && mv tmp noid_genes.gff


### First I format SweeD outputs and add the chrom/scaffold name to all files, while removing first 3 lines (empty line, garbage and header)
########## CHANGE PREFIX OF SWEED REPORT FILE FOR EACH SPECIES BEFORE RUNNING THIS SCRIPT ##########

while read p; do

tail -n +3 OmegaPlus_Report.weigel_athaliana_IBE.vcf.gz_$p| awk -v var="$p" 'BEGIN {OFS="\t"} {print var, int($1), int($1), "CLRSCAN", $2}' > formatted\_$p

done < chrom.txt




### Concatenate formatted outputs into 1 and remove all the formatted chrom/scaffolds

cat formatted_* > all.bed
rm formatted*


### Restrict analysis to genes only

/data/programs/bedtools2/bin/bedtools intersect -a all.bed -b noid_genes.gff -wa | uniq > tmp
mv tmp all.bed


### Here I make a list of genes intersecting the clr scans. Each gene gets reported as many times as the number of clr scan it includes

/data/programs/bedtools2/bin/bedtools intersect -a noid_genes.gff -b all.bed -wo | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$14}' > genes_stats.gff



### Here first I add for each gene, the count of CLR scans - then I format the final file to only have gene name (given by chrom:start-end), CLR and number of CLR scans within the gene
### I also remove the 500bp flanks from each side of genes, as I need unflanked genes to match genes with orthogroups using the map I created from Jim file

/data/programs/bedtools2/bin/bedtools coverage -a genes_stats.gff -b all.bed -counts | awk '{OFS="\t"}{if ($4 > 1) print $1":"$4+500"-"$5-500,$9,$10; else print $1":"$2"-"$5-500,$9,$10;}' > final_genes_analysis.txt


Rscript results_analysis.R

rm final_genes_analysis.txt
rm genes_stats.gff
rm all*
rm noid_genes.gff
