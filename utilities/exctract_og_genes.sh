for file in ./OG*
do
	sed 's/:/___/g' $file > tmp.txt

	while read line
		do
    		grep -A 1 -h $line *faa >> $file\_extracted_genes.txt
    		# add commands to process the line here
	done < tmp.txt


done


rm tmp.txt
