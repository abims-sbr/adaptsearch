#/bin/bash

for file in *.fasta 
do      
	sed -i -e 's/^.*[ ]([+])[ ]/>/g' -e 's/[:].*$//' $file
	#sed -i -E 's/^>.+\(\+\) ([^:]+):.+$/>\1/' $file #V2
done
