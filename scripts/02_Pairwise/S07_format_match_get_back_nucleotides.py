#!/usr/bin/env python
## AUTHOR: Eric Fontanillas
## LAST VERSION: 14/08/14 by Julie BAFFARD

## DESCRIPTION: find back the nucleotide sequence corresponding to the protein sequences in file "09_onlyMatch_filtered_200bp.fasta"

## 1 ARGUMENT:  Minimum length (e.g. 100)
#MIN_LENGTH = 100

## SPLIT for string.split fasta name
SPLIT = "||"


############################
## DEF 1 : Generates bash ## 
############################
## key = fasta name; value = sequence (WITH GAP, IF ANY, REMOVED IN THIS FUNCTION)
def dico(fasta_file):
    count_fastaName=0
    bash1 = {}
    with open(fasta_file, "r") as F1:        
        for line, line2 in itertools.izip_longest(*[F1]*2):            
            if line[0] == ">":
                count_fastaName = count_fastaName + 1
                fasta_name = line[1:-1]                
                sequence = line2[:-1]
                bash1[fasta_name] = sequence
    return(bash1, count_fastaName)
#####################################


#####################################
## DEF 2 : Get list of fasta names ##
#####################################
def get_list_fasta_names(bash):
    list_fasta_name = []
    L = bash.keys()

    for name in L:
        S1 = string.split(name, SPLIT) 
        short_name = S1[0]
        long_name = name
        L = [short_name, long_name]
        
        list_fasta_name.append(L)

    return(list_fasta_name)
#################################


#####################################################
## DEF 3 :Get sequences from a list of fasta_names ##
#####################################################
def get_sequences(list_fasta_names, bash2):
    L = bash2.keys()
    ln = len(list_fasta_names)
    bash3 = {}

    i = 1
    for names in list_fasta_names:
        short_name = names[0]
        long_name = names[1]
        
        for fasta_name in L:
            if short_name in fasta_name:
                new_fasta_name = long_name
                new_fasta_seq = bash2[fasta_name]
                bash3[new_fasta_name] = new_fasta_seq
        i = i+1
    return(bash3)
#####################################


#####################
#### RUN RUN RUN ####
#####################
import string, os, sys, itertools

## 1 ## INPUT/OUTPUT
SHORT_FILE = sys.argv[2] #short-name-query_short-name-db

path_IN1 = "%s/09_onlyMatch_filtered_%s.fasta" %(SHORT_FILE, SHORT_FILE) 
path_IN2 = sys.argv[1]  ## Initially the DB in blast input

path_OUT = "%s/09_onlyMatch_filtered_nucleotideBACK_%s.fasta" %(SHORT_FILE, SHORT_FILE)
file_OUT = open(path_OUT, "w")

## 2 ## RUN
bash1, count_FastaName1 = dico(path_IN1)            ### DEF1 ###
ln1 = len(bash1.keys())

bash2, count_FastaName2  = dico(path_IN2)           ### DEF1 ###
ln2 = len(bash2.keys())

list_fasta_names = get_list_fasta_names(bash1)      ### DEF2 ###
ln = len(list_fasta_names)

bash3 = get_sequences(list_fasta_names, bash2)      ### DEF3 ###
for fasta_name in bash3:    
    file_OUT.write(">%s\n" %fasta_name)
    file_OUT.write("%s\n" %bash3[fasta_name])    


file_OUT.close()
