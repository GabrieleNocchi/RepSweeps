#!/bin/bash
for file in *.txt; do
  sed -i '1,6d' "$file"
  sed -i 's/-nan/0/g' "$file"
  result=$(awk 'NR>1{for(i=2;i<=NF-2;i++) if($i!=0) {sum += $i; count++}} END{print sum/count}' "$file")
  echo "$file $result" >> output.txt
done
