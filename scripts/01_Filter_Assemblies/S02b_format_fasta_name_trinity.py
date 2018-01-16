#!/usr/bin/env python
## AUTHOR: Eric Fontanillas
## LAST VERSION: 06.12.2011
## DESCRIPTION: format fasta name in TRINITY output

from os import listdir
import re

###################
###### DEF 1 ######
###################
def dico_format_fasta_name(path_in, SUFFIX):
    f_in = open(path_in, "r")
    bash = {}
    file_read = f_in.read()
    S1 = file_read.split(">")
    k = 0

    for element in S1:
        if element != "":
            S2 = element.split("\n")
            fasta_name = S2[0]
            fasta_seq = S2[1]
            L = fasta_name.split("_")
            match=re.search('(\D+)(\d+)', L[0])
            short_fasta_name= SUFFIX + match.group(2) + "_" + L[1] + "_" + L[2]
            bash[short_fasta_name] = fasta_seq

    return bash
#~#~#~#~#~#~#~#~#~#

###################
### RUN RUN RUN ###
###################
import string, os, sys, re

path_IN = sys.argv[1]
path_OUT = sys.argv[2]
suffix= sys.argv[3]
file_OUT = open(path_OUT, "w")
#Extract suffix info

dico = dico_format_fasta_name(path_IN, suffix)   ### DEF1 ###

print((len(list(dico.keys()))))

KB = list(dico.keys())

## Sort the fasta_name depending their number XX : ApXX
BASH_KB = {}
for name in KB:    
    L = name.split("_")
    nb = L[0][2:]    
    nb = int(nb)    
    BASH_KB[nb] = name

KKB = list(BASH_KB.keys())
KKB.sort()

for nb in KKB:
    fasta_name = BASH_KB[nb]
    seq = dico[fasta_name]
    file_OUT.write(">%s\n" %fasta_name)
    file_OUT.write("%s\n" %seq)

file_OUT.close()
