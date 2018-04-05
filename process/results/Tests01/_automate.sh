#!/usr/bin/bash

LIST=$(ls -1 *.csv)

echo $LIST

for _name in $LIST; do
	awk 'BEGIN{CONVFMT="%.9f"; FS=OFS=" "}{for(i=1; i<=NF; i++)if($i~/^[0-9]+([eE][+-][0-9]+)?/)$i+=0;}1' $_name > temp.txt
	cat temp.txt > $_name
done


