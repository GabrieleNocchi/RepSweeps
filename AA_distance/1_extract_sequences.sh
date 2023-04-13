for file in ./*.txt
do
	sed 's/:/___/g' $file > tmp

	while read line
		do
    		
	grep -A 1 -h $line /data/home/gabriele/RepSweeps/SweeD_run_I/aa_distance/*faa >> $file\_extracted.fasta
    		# add commands to process the line here
	done < tmp
	
done


rm tmp
