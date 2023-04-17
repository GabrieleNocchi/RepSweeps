#!/bin/bash
for file in *.txt; do
  sed -i '1,6d' "$file"
  result=$(awk 'NR>1{for(i=2;i<=NF-2;i++) if($i!="nan") {sum += $i; count++}} END{print sum/count}' "$file")
  echo "$file $result" >> output.txt
done

