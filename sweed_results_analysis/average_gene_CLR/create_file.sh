while read p; do
        grep -P "$p" ortho_gene_count.txt
done < ortho_list.txt
