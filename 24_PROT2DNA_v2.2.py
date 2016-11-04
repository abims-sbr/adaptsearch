#!/usr/local/public/python-2.6.5/bin/python
#/usr/bin/python ==> PYTHON 2.4

## Author: Eric FONTANILLAS
## Last modification: 12/11/10
## Object:


############################
##### DEF1 : Get Pairs #####
############################
def get_pairs(fasta_file_path):
    F2 = open(fasta_file_path, "r")
    list_pairwises = []
    #bash_pairwises = {}
    #bash_pairwises_big = {}
    
    h = 0
    while 1:
        next2 = F2.readline()
        if not next2:
            break
        if next2[0] == ">":
            ## Read 1 pairwise
            fasta_name_query = next2[1:-1]
            next3 = F2.readline()
            fasta_seq_query = next3[:-1]
            next3 = F2.readline()    ## jump one empty line (if any after the sequence)
            fasta_name_match = next3[1:-1]
            next3 = F2.readline()
            fasta_seq_match = next3[:-1]
            h = h+1

            ## Format shorter names 
            S1 = string.split(fasta_name_query, "||")
            short_fasta_name_query = S1[0]

            S2 = string.split(fasta_name_match, "||")
            short_fasta_name_match = S2[0]
           
            ## create and record pairwise (in list format)
            pairwise = [short_fasta_name_query,fasta_seq_query ,short_fasta_name_match, fasta_seq_match]
            list_pairwises.append(pairwise)

    F2.close()

    print "Number of pairwises parsed = %d" %h
    return(list_pairwises)
##############################################

#####################################################################
##### DEF1_modified : Get Pairs (only fasta names, not squences #####
#####################################################################
def get_pairs_modified(fasta_file_path):
    F2 = open(fasta_file_path, "r")
    list_pairwises = []
    list_query = []
    list_match= []
    #bash_pairwises = {}
    #bash_pairwises_big = {}
    
    h = 0
    while 1:
        next2 = F2.readline()
        if not next2:
            break
        if next2[0] == ">":
            ## Read 1 pairwise
            fasta_name_query = next2[1:-1]
            next3 = F2.readline()
            fasta_seq_query = next3[:-1]
            next3 = F2.readline()    ## jump one empty line (if any after the sequence)
            fasta_name_match = next3[1:-1]
            next3 = F2.readline()
            fasta_seq_match = next3[:-1]
            h = h+1

            ## Format shorter names 
            S1 = string.split(fasta_name_query, "||")
            short_fasta_name_query = S1[0]

            S2 = string.split(fasta_name_match, "||")
            short_fasta_name_match = S2[0]
           
            ## create and record pairwise (in list format)
            pairwise = [short_fasta_name_query,short_fasta_name_match]
            list_pairwises.append(pairwise)
            list_query.append(short_fasta_name_query)
            list_match.append(short_fasta_name_match)
            

    F2.close()

    print "Number of pairwises parsed = %d" %h
    return(list_pairwises, list_query, list_match)
##############################################

##################
###### DEF2 ######
##################
## Generates bash, with key = fasta name; value = sequence (WITH GAP, IF ANY, REMOVED IN THIS FUNCTION)

def dico(fasta_file):

    count_fastaName=0
    F1 = open(fasta_file, "r")
    
    bash1 = {}
    while 1:
        nextline = F1.readline()
        
        if not nextline :
            break
        
        if nextline[0] == ">":
            count_fastaName = count_fastaName + 1
            fasta_name = nextline[1:-1]
            nextline = F1.readline()
            sequence = nextline[:-1]
            
            if fasta_name not in bash1.keys():
                bash1[fasta_name] = sequence
            else:
                print fasta_name
                

    F1.close()
    return(bash1)
#####################################


##################
###### DEF3 ######
##################

def get_dico_seq_subset(path, list_fastaNames):
    F1 = open(path, "r")

    nb_line_treated = 0
    
    bash_subset = {}
    
    while 1:
        nextline = F1.readline()

        nb_line_treated = nb_line_treated + 1
        if nb_line_treated%10000 == 0:
            print "\t%d" %nb_line_treated
        
        if not nextline :
            break
        
        if nextline[0] == ">":
            if nextline[1:-1] in list_fastaNames:
                name = nextline[1:-1]
                sequence = F1.readline()
                sequence = sequence[:-1]
                bash_subset[name] = sequence
    return(bash_subset)
#####################################

########################
##### RUN RUN RUN ######
########################

import sys, string, os, re, copy, itertools

## 1 ## INPUT/OUTPUT
path1 = sys.argv[1]  ## initial Query in the RBH
path2 = sys.argv[2]  ## initial DB in the RBH
WORK_DIR = sys.argv[3]

path3 = "%s/19_ReciprocalBestHits.fasta" %WORK_DIR

file_OUT = open("%s/25_DNAalignment_corresponding_to_protein_from_19_RBH.fasta" %WORK_DIR, "w")

## 2 ## Get DB
print "Get list of fasta name involved in RBH"
list_pairwise,list_query,list_match = get_pairs_modified(path3) ### DEF1_modified ###

print "Get subset of Alvinella db"
bash_subset1 = get_dico_seq_subset(path1, list_query)   ### DEF3 ###

print "Get subset of Paralvinella db"
bash_subset2 = get_dico_seq_subset(path2, list_match)   ### DEF3 ###


## 3 ## Grab Pairwise DNA
for pair in  list_pairwise:  # list pairwise comes from '19_RBH ...'
    name1 = pair[0]     # Ap
    seq1 = bash_subset1[name1]
    #print seq1
    name2 = pair[1]     # Ps
    seq2 = bash_subset2[name2]
    #print seq2

    file_OUT.write(">%s\n" %name1)
    file_OUT.write("%s\n" %seq1)
    file_OUT.write(">%s\n" %name2)
    file_OUT.write("%s\n" %seq2)
    
file_OUT.close()
