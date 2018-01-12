#/bin/bash

# v2 - this script modifies the 'Orthogroups.txt' file in order to make it easily readable by the following script, filter_orthofinder.py
  #Example :
    #OG0000001: Gene.117__As119_1/1_1.000_543__g.117__m.117 Gene.157__As170_1/1_1.000_1203__g.157__m.157
  #Becomes :
    #As119_1/1_1.000_543 As170_1/1_1.000_1203
    
# removes 'OGxxxxxxx: '
sed -E 's/OG[0-9]{7,}:\s//' $1 > $2
# removes things like Gene.119__
sed -i -E 's/Gene\.[0-9]{1,}\_\_/>/g' $2
# removes things like __g.117__m.117
sed -i -E 's/\_\_g\.[0-9]{1,}\_\_m\.[0-9]{1,}//g' $2

# Old version

# removes 'OGxxxxxxx '
#sed -E 's/OG[0-9]{7}:\s//' $1 > $2
# replace _+_ by (+) because '_' causes bugs
#sed -i 's/_+_/(+)/g' $2
# Replaces everything by '>'
#sed -i -E 's/m\.[0-9]{1,}[^()]+\(\+\)\s*/>/g' $2
# Removes terminal '(+)'
#sed -i 's/(+)//g' $2
# Removes last suite of unwanted numbers, underscore and dash
#sed -i -E 's/\_[0-9]{1,}-[0-9]{1,}//g' $2
