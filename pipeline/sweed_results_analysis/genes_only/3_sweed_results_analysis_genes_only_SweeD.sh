regex='(.+)	(.+)'


while read p; do

if [[ $p =~ $regex ]];  then tail -n +4 SweeD_Report.kubota_ahalleri.vcf.gz_${BASH_REMATCH[2]} | awk -v var="${BASH_REMATCH[2]}" 'BEGIN {OFS="\t"} {print var, $2, 5}' > formatted_${BASH_REMATCH[2]}; fi

done < genes_coord.txt



cat formatted_* > final_genes_analysis.txt
rm formatted*


Rscript results_analysis.R
