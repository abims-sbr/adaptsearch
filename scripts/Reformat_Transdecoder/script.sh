#/bin/bash

#This script contains regex to re-write the outputs of transdecoder to the original AdaptSearch format 

# removes 'Gene.xxxx::'
sed -E 's/Gene.[0-9]+:://g' $1 > $2
# remove eveything from the first ':' until the end of the line ('$')
sed -i -E 's/:.*$//g' $2

python oneLine.py $2 "Formatted_"$1