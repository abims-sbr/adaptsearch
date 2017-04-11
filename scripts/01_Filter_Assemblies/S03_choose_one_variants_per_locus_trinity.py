#!/usr/bin/env python

## AUTHOR: Eric Fontanillas

## LAST VERSION: 06.12.2011

## DESCRIPTION: Remove redondant transcripts (i.e. transcript from the same locus) from TRINITY on the basis of 1 criteria: 
            ## 1. [CRITERIA 1] choose the longuest sequence (once any "N" have been removed => effective_length = length - N number

###################
###### DEF 1 ######
###################
def dico_filtering_redundancy(path_in):
    f_in = open(path_in, "r")
    bash = {}
    bash_unredundant = {}
    file_read = f_in.read()
    S1 = string.split(file_read, ">")
    k = 0


    ## 1 ## Extract each transcript and group them in same locus if they share the same "short_fasta_name"
    for element in S1:
        if element != "":
            S2 = string.split(element, "\n")
            fasta_name = S2[0]
            fasta_seq = S2[1]

            L = string.split(fasta_name, "_")

            ## 1.1. ## Extract short fasta name
            short_fasta_name = L[0] + L[1]

            #new_fasta_name = L[0] + "_" + L[1] + "_" + L[2]
            
            ## 1.2. ## Extract length
            # Length extraction from the name not necessary
            #lgt = L[-2]
            #list_lgt = string.split(lgt, ":")
            #length = string.atoi(list_lgt[1])
            #print "\tLength in fasta name = %.5f" %length
            
            ## Used later for [CRITERIA 1] (see below)
            countN = string.count(fasta_seq, "N")
            #print "\tNb of N in sequence = %d" %countN
            length = len(fasta_seq)
            effective_length = length - countN
            #####################################################

            if short_fasta_name not in bash.keys():
                bash[short_fasta_name] = [[fasta_name, fasta_seq, effective_length]]
            else:
                bash[short_fasta_name].append([fasta_name, fasta_seq, effective_length])
        k = k+1
        if k%1000 == 0:
            print (k)
    f_in.close()

    for key in bash.keys():
        #print "%s = %d" %(key,len(bash[key]))

        ## 2 ## IF ONE TRANSCRIPT PER LOCUS:
        ## In this case => we record directly
        if len(bash[key]) == 1:
            #print bash[key]
            entry = bash[key][0]
            name = entry[0]
            seq = entry[1]
            bash_unredundant[name] = seq

        ## 3 ## IF MORE THAN ONE TRANSCRIPTS PER LOCUS:
        ## In this case:
        ## [CRITERIA 1]: Choose the longuest sequence (once any "N" have been removed => effective_length = length - N numb
        elif len(bash[key]) > 1:   ### means there are more than 1 seq
            MAX_LENGTH = {}
            for entry in bash[key]:    ## KEY = short fasta name    || VALUE = list of list, e.g. :  [[fasta_name1, fasta_seq1],[fasta_name2, fasta_seq2][fasta_name3, fasta_seq3]]
                #print entry
                name = entry[0]
                seq = entry[1]
                effective_length = entry[2]

                ## Bash for [CRITERIA 1]
                MAX_LENGTH[effective_length] = entry

            
            ## Sort keys() for MAX_LENGTH bash 
            KC = MAX_LENGTH.keys()
            KC.sort()
            
            ## Select the best entry
            MAX_LENGTH_KEY = KC[-1]  ## [CRITERIA 1]
            BEST_ENTRY = MAX_LENGTH[MAX_LENGTH_KEY]
                
            BEST_fasta_name = BEST_ENTRY[0]
            BEST_seq = BEST_ENTRY[1]
            bash_unredundant[BEST_fasta_name] = BEST_seq
        
    return bash_unredundant      
#~#~#~#~#~#~#~#~#~#



###################
### RUN RUN RUN ###
###################
import string, os, sys, re

path_IN = sys.argv[1]
path_OUT = sys.argv[2]
file_OUT = open(path_OUT, "w")
#length_seq_max=sys.argv[3]
dico = dico_filtering_redundancy(path_IN)   ### DEF1 ###

#print dico.keys()

#print len(dico.keys())

KB = dico.keys()
#print "dis donc %s" %KB

## Sort the fasta_name depending their number XX : ApXX
BASH_KB = {}
for name in KB:
    #print name
    L = string.split(name, "_")
    nb = L[0][2:]
    nb = string.atoi(nb)
    #print nb
    BASH_KB[nb] = name

KKB = BASH_KB.keys()
KKB.sort()

for nb in KKB:
    fasta_name = BASH_KB[nb]
    #print fasta_name
    seq = dico[fasta_name]
    #if int(len(seq)) > int(length_seq_max):
    file_OUT.write(">%s\n" %fasta_name)
    file_OUT.write("%s\n" %seq)

file_OUT.close()
