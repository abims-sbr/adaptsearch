#/bin/bash


#This script contains regex to re-write the outputs of transdecoder to the original AdaptSearch format 
#Example :
#OG0007971: m.35 g.35  ORF g.35 m.35 type_internal len_307 _+_ Th132_1/1_1.000_923_1-924_+_
#Becomes :
#Th132_1/1_1.000_923

# removes 'OGxxxxxxx '
sed -i -E 's/OG[0-9]{7}:\s//' $1 
# replace _+_ by (+) because '_' causes bugs
sed -i 's/_+_/(+)/g' $1
# Replaces everything by '>'
sed -i -E 's/m\.[0-9]{1,}[^()]+\(\+\)\s*/>/g' $1
# Removes terminal '(+)'
sed -i 's/(+)//g' $1
# Removes last suite of unwanted numbers, underscore and dash
sed -i -E 's/\_[0-9]{1,}-[0-9]{1,}//g' $1
