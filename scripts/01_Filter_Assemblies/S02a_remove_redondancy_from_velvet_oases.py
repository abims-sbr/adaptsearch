#!/usr/bin/env python

## AUTHOR: Eric Fontanillas

## LAST VERSION: 09.06.2011

## DESCRIPTION: Remove redondant transcripts (i.e. transcript from the same locus) from Oases output on the basis of two recursive criterias (see in DEF1):  
            ## 1. [CRITERIA 1] Keep in priority seq with BEST "confidence_oases_criteria" present in the fasta name
            ## 2. [CRITERIA 2] Second choice (if same coverage) : choose the longuest sequence (once any "N" have been removed => effective_length = length - N number
## => criticize of this approach: the transcripts may come from a same locus but may be not redundant (non-overlapping) ==> SEE "DEF2" for an alternative

###################
###### DEF 1 ######
###################
def dico_filtering_redundancy(path_in):
    f_in = open(path_in, "r")
    bash = {}
    bash_unredundant = {}
    file_read = f_in.read()
    # S1 = string.split(file_read, ">")
    S1 = file_read.split(">")
    k = 0


    ## 1 ## Extract each transcript and group them in same locus if they share the same "short_fasta_name"
    for element in S1:
        if element != "":
            # S2 = string.split(element, "\n")
            S2 = element.split("\n")
            fasta_name = S2[0]
            fasta_seq = S2[1:-1] # that line was unindented
            fasta_seq = "".join(fasta_seq) # that line was unindented
            # L = string.split(fasta_name, "_")
            L = fasta_name.split("_")
            short_fasta_name = L[0] + L[1]
            
            
            #####################################################
            ## Used later for [CRITERIA 1] (see below)
            confidence_oases_criteria = L[-3]

            ## Used later for [CRITERIA 1] (see below)
            # countN = string.count(fasta_seq, "N")
            countN = fasta_seq.count("N")
            length = len(fasta_seq)
            effective_length = length - countN
            #####################################################
            if short_fasta_name not in list(bash.keys()): 
                bash[short_fasta_name] = [[fasta_name, fasta_seq, confidence_oases_criteria, effective_length]]
            else:
                bash[short_fasta_name].append([fasta_name, fasta_seq, confidence_oases_criteria, effective_length])
        k = k+1
        f_in.close()

    for key in list(bash.keys()):
        ## 2 ## IF ONE TRANSCRIPT PER LOCUS:
        ## In this case => we record directly
        if len(bash[key]) == 1:
            entry = bash[key][0]
            name = entry[0]
            seq = entry[1]
            bash_unredundant[name] = seq

        ## 3 ## IF MORE THAN ONE TRANSCRIPTS PER LOCUS:
        ## In this case:
        ## 1. [CRITERIA 1] Keep in priority seq with BEST "confidence_oases_criteria" present in the fasta name
        ## 2. [CRITERIA 2] Second choice (if same coverage) : choose the longuest sequence (once any "N" have been removed => effective_length = length - N numb
        elif len(bash[key]) > 1:   ### means there are more than 1 seq            
            MAX_CONFIDENCE = {}
            MAX_LENGTH = {}
            for entry in bash[key]:    ## KEY = short fasta name    || VALUE = list of list, e.g. :  [[fasta_name1, fasta_seq1],[fasta_name2, fasta_seq2][fasta_name3, fasta_seq3]]
                name = entry[0]
                seq = entry[1]
                effective_length = entry[3]
                confidence_oases_criteria = entry[2]

                ## Bash for [CRITERIA 2]
                MAX_LENGTH[effective_length] = entry

                ## Bash for [CRITERIA 1]
                # confidence_oases_criteria = string.atof(confidence_oases_criteria)
                confidence_oases_criteria = float(confidence_oases_criteria)
                if confidence_oases_criteria not in list(MAX_CONFIDENCE.keys()):
                    MAX_CONFIDENCE[confidence_oases_criteria] = entry
                else:    ## IF SEVERAL SEQUENCES WITH THE SAME CONFIDENCE INTERVAL => RECORD ONLY THE LONGUEST ONE [CRITERIA 2]
                    current_seq_length = effective_length
                    yet_recorded_seq_length = MAX_CONFIDENCE[confidence_oases_criteria][3]
                    if current_seq_length > yet_recorded_seq_length:
                        MAX_CONFIDENCE[confidence_oases_criteria] = entry   ## Replace the previous recorded entry with the same confidence interval but lower length
            
            
            
            ## Sort keys() for MAX_CONFIDENCE bash 
            KC = list(MAX_CONFIDENCE.keys())
            KC.sort()
            
            ## Select the best entry
            MAX_CONFIDENCE_KEY = KC[-1]  ## [CRITERIA 1]
            BEST_ENTRY = MAX_CONFIDENCE[MAX_CONFIDENCE_KEY]
            # if len(KC) > 1:    ## [CRITERIA 1] Different confidence criteria found
            #     MAX_CONFIDENCE_KEY = KC[-1]
            #     BEST_ENTRY = MAX_CONFIDENCE[MAX_CONFIDENCE_KEY]
            # else:             ## [CRITERIA 2] ALL transcripts have the same confidence interval => we should use the second criteria: effective length (= length - nb of N)
            #     MAX_LENGTH_KEY = KL[-1]
            #     BEST_ENTRY = MAX_LENGTH[MAX_LENGTH_KEY] 
                
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
#length_seq_max=sys.argv[3]
file_OUT = open(path_OUT, "w")

dico = dico_filtering_redundancy(path_IN)   ### DEF1 ###


KB = list(dico.keys())

## Sort the fasta_name depending their number XX : ApXX
BASH_KB = {}
for name in KB:
    # L = string.split(name, "_")
    L = name.split("_")
    # nb = string.atoi(L[1])
    nb = int(L[1])
    BASH_KB[nb] = name
NEW_KB = []    
KKB = list(BASH_KB.keys())
KKB.sort()

for nb in KKB:
    fasta_name = BASH_KB[nb]
    seq = dico[fasta_name]
    #if int(len(seq)) > int(length_seq_max):
    file_OUT.write(">%s\n" %fasta_name)
    file_OUT.write("%s\n" %seq)

file_OUT.close()
