#!/usr/bin/python

## AUTHOR: Eric Fontanillas
## LAST VERSION: 20/08/14 by Julie BAFFARD

## DESCRIPTION: formatting a fasta format into phylip format for using with PAML

import string, os, sys

if len(sys.argv) == 1:
    print "put arguments!!"
    print "USAGE: $phylip2fasta.py INPUT OUTPUT"


## INPUT
f1 = sys.argv[1]
F1 = open("%s" %f1, 'r')

## OUTPUT
f2 = sys.argv[2]
F2 = open("%s" %f2, 'w')

###### def1 ######
# Dans un multialignement fasta, cette fonction permet de formatter les noms de chaque sequence fasta

def format(File_IN):
    c = 0
    fichier = ""
    while 1 :
        c = c + 1
        next = File_IN.readline()
        if not next :
            break
        
        S1 = string.split(next, " ")
        fasta_name = S1[0]
        fasta_seq = S1[1][:-1]
        fichier = fichier + ">" + fasta_name + "\n" + fasta_seq + "\n"
        
    return (fichier,c)
#-#-#-#-#-#-#-#-#-#-#

###################
### RUN RUN RUN ###
###################

F1.readline() ## jump the first line

fichier_txt, c = format(F1)   ### DEF1 ###

F2.write(fichier_txt)

F1.close()
F2.close()
