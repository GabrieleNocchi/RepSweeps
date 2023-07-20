Rscript extract_all_genes.R

awk '{print $1}' all_genes_average.txt | sed 's/:\|-/\t/g' | tail -n +2 > genes.bed


tail -n +2 hannus_recomb.txt | awk '{OFS="\t"}{print $1,$2,$3,$6}' > edited_recomb.txt






/data/programs/bedtools2/bin/bedtools intersect -a genes.bed -b edited_recomb.txt -wo | awk '{OFS="\t"}{print $1,$2,$3,$7}' > genes_all_rec.bed


/data/programs/bedtools2/bin/bedtools coverage -a genes_all_rec.bed -b edited_recomb.txt -counts | awk '{OFS="\t"}{print $1":"$2"-"$3,$4,$5}' > final_genes_all_rec.txt


Rscript make_rec_emp.R


rm final_genes_all_rec.txt
rm genes_all_rec.bed
rm edited_recomb.txt
rm genes.bed
