regex='(.+)	(.+)'


while read p; do

if [[ $p =~ $regex ]];  then tail -n +3 OmegaPlus_Report.kubota_ahalleri.vcf.gz_${BASH_REMATCH[2]} | awk -v var="${BASH_REMATCH[2]}" 'BEGIN {OFS="\t"} {print var, $2, 5}' > formatted_${BASH_REMATCH[2]}; fi

done < genes_coord.txt



cat formatted_* > final_genes_analysis.txt
rm formatted*
awk 'NR % 5 == 3' final_genes_analysis.txt > final_genes_analysis_1.txt

Rscript results_analysis.R
