awk '{OFS="\t"}NR>1{print $1,$4-1000,$5+1000}' genes.gff | awk '{OFS="\t"}{if ($2 > 0)  print $1":"$2"-"$3; else print $1":"1"-"$3;}' > genes_coord.txt
tail -n +2 genes.gff | awk '{print $1":"$4"-"$5}' > g
paste genes_coord.txt g > merged
rm g
mv merged genes_coord.txt
split --numeric-suffixes=1 -n l/999 -d genes_coord.txt splitted
