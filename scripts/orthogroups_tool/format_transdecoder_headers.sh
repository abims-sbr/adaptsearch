#/bin/bash

# This script contains regex to re-write the outputs of transdecoder to the original AdaptSearch format 

sed -i -E 's/OG[0-9]{7}:\s//' $1
sed -i 's/_+_/(+)/g' $1
sed -i -E 's/m\.[0-9]{1,}[^()]+\(\+\)\s*/>/g' $1
#; s/:\S+//g
sed -i -E 's/_[0-9]{1,}-[0-9]{1,}\(\+\)//g' $1
