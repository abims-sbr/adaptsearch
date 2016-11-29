#!/usr/bin/python

## AUTHOR: Eric Fontanillas

## LAST VERSION: 06.12.2011

## DESCRIPTION: format fasta name in TRINITY output

SUFFIX = "Pp"


###################
###### DEF 1 ######
###################
def dico_format_fasta_name(path_in, SUFFIX):
    f_in = open(path_in, "r")
    bash = {}
    file_read = f_in.read()
    S1 = string.split(file_read, ">")
    k = 0

    for element in S1:
        if element != "":
            S2 = string.split(element, "\n")
            fasta_name = S2[0]
            fasta_seq = S2[1]

            L = string.split(fasta_name, "_")

            short_fasta_name = L[0] + "_" + L[1] + "_" + L[2]

            short_fasta_name = string.replace(short_fasta_name, "comp", SUFFIX)
            
            bash[short_fasta_name] = fasta_seq
        
    return bash
#~#~#~#~#~#~#~#~#~#



###################
### RUN RUN RUN ###
###################
import string, os, sys, re

path_IN = sys.argv[1]
path_OUT = sys.argv[2]
file_OUT = open(path_OUT, "w")

dico = dico_format_fasta_name(path_IN, SUFFIX)   ### DEF1 ###


print len(dico.keys())

KB = dico.keys()

## Sort the fasta_name depending their number XX : ApXX
BASH_KB = {}
for name in KB:
    print name
    
    L = string.split(name, "_")
    nb = L[0][2:]
    nb = string.atoi(nb)
    print nb
    BASH_KB[nb] = name

KKB = BASH_KB.keys()
KKB.sort()

for nb in KKB:
    fasta_name = BASH_KB[nb]
    seq = dico[fasta_name]
    file_OUT.write(">%s\n" %fasta_name)
    file_OUT.write("%s\n" %seq)

file_OUT.close()
