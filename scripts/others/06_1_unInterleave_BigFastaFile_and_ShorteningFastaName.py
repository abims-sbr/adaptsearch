#!/usr/bin/python

## AUTHOR: Eric Fontanillas

## LAST VERSION: 18.05.2011

## DESCRIPTION: Uninterleaved transcript file ("transcripts.fa", i.e. the oases output)  => create the file "transcripts.fasta"  and format the fasta name (to make it shorter)



###############################
### DEF : FORMAT FASTA NAME ###
###############################
def format_fastaName(fastaName):
    # INPUT: EXAMPLE for fastaNAME = Locus_25599_Transcript_1/1_Confidence_1.000_Length_215
    # OUTPUT: We want: Ac25599_1/1_1.000_215
    l = string.split(fastaName, "_")
    fastaName = SPECIES + l[1] + "_" + l[3]+ "_" + l[5]+ "_" + l[7] 
                                         
    return fastaName

###########
## DEF 1 ##
###########
## Generates bash, with key = fasta name; value = sequence (WITH GAP, IF ANY, REMOVED IN THIS FUNCTION)

def uninterleave_BigFastaFile(path_fileIN, path_fileOUT):
    F1 = open(path_fileIN, "r")
    F2 = open(path_fileOUT, "w")

    bash1 = {}
    j = 0
    seq = ""
    k = 0
    while 1:
        k = k+1
        if k%10000 ==0:
            print "\t%d" %k
            
        nextline = F1.readline()
        
        if not nextline :
            seq = string.replace(seq, "\n", "")
            ## record the last entrie
            F2.write(">%s\n" %fasta_name)
            F2.write("%s\n" %seq)
            break
        
        if nextline[0] != ">":
            seq = seq+nextline
            j = 1

        elif nextline[0] == ">":
            ## 1 ## record previous sequence
            if j ==1:
                seq = string.replace(seq, "\n", "")
                F2.write(">%s\n" %fasta_name)
                F2.write("%s\n" %seq)
                
            ## 2 ## new fasta name
            seq = ""            
            fasta_name = nextline[1:-1]
            fasta_name = format_fastaName(fasta_name)
            print fasta_name
            
            
    F1.close()
    F2.close()
    return()
#####################################

import string, os, sys

#SPECIES = "Ac_"
SPECIES = sys.argv[1]    ## format "Ac_"


pathIN = "../tmp"
List_onlydirectories = [name for name in os.listdir(pathIN) if os.path.isdir(os.path.join(pathIN, name))]
List_onlydirectories.sort()

for kmerRun in List_onlydirectories:
    print kmerRun

    path_fileIN = "%s/%s/transcripts.fa" %(pathIN, kmerRun)
    path_fileOUT = "%s/%s/transcripts.fasta" %(pathIN, kmerRun)
    
    uninterleave_BigFastaFile(path_fileIN, path_fileOUT)    ### DEF1 ###
    
    
