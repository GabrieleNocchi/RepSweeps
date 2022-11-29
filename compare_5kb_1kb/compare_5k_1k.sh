while read p; do

tail -n +4 SweeD_Report.savolainen_alyrata_scandinavia.vcf.gz_$p| awk -v var="$p" 'BEGIN {OFS="\t"} {print var, $1-5000, $1+5000, "CLRSCAN", $2}' | awk '{OFS="\t"}{if ($2 > 0)  print $1, $2, $3, $4, $5;else print $1, 0, $3, $4, $5;}' > formatted\_$p

done < chrom.txt




### Concatenate formatted outputs into 1 and remove all the formatted chrom/scaffolds
cat formatted_* > all_5k.bed

rm formatted*



/data/programs/bedtools2/bin/bedtools coverage -a all_5k.bed -b savolainen_alyrata_scandinavia.vcf.gz -counts > tmp
mv tmp all_5k.bed

/data/programs/bedtools2/bin/bedtools intersect -a all_5k.bed -b all_1k.bed -wo | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$11}' | awk '{OFS="\t"}{print $1":"$2"-"$3,$2,$3,$4,$5,$6,$7}'> to_plot_in_R.txt



