#!/usr/local/public/python-2.6.5/bin/python
#/usr/bin/python ==> PYTHON 2.4

## AUTHOR: Eric Fontanillas
## LAST VERSION: 14/08/14 by Julie BAFFARD


#####################################################################
##### DEF1_modified : Get Pairs(only fasta names, not squences) #####
#####################################################################
def get_pairs_modified(fasta_file_path):
    F2 = open(fasta_file_path, "r")
    list_pairwises = []
    list_query = []
    list_match= []
    
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

    print "Number of pairwises parsed = %d \n\n\n" %h
    return(list_pairwises, list_query, list_match)
##############################################


##################
###### DEF2 ######
##################
def get_dico_seq_subset(path, list_fastaNames):
    F1 = open(path, "r")

    nb_line_treated = 0
    
    bash_subset = {}
    
    while 1:
        nextline = F1.readline()

        nb_line_treated = nb_line_treated + 1
        
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
SHORT_FILE = sys.argv[3]  ## short-name-query_short-name-db

### for print
name1 = sys.argv[4]
name2 = sys.argv[5]
option_e = string.atof(sys.argv[6])

path3 = "%s/19_ReciprocalBestHits_%s.fasta" %(SHORT_FILE, SHORT_FILE)
os.system("cp -fr %s ./" %path3)

file_OUT = open("./25_DNAalignment_corresponding_to_protein_from_19_RBH_%s.fasta" %SHORT_FILE, "w")

## 2 ## RUN
##Get DB
##Print
print "-------------------- Pairwise %s --------------------\n" %SHORT_FILE
print "database : %s" %name1
print "query file : %s" %name2
print "\n***** RUN FIRST BLAST *****\n\n"
print "database : %s" %name2
print "query file : only the sequences of %s which matched during the last BLAST" %name1
print "\n***** RUN SECOND BLAST *****\n\n"

#Get list of fasta name involved in RBH
list_pairwise,list_query,list_match = get_pairs_modified(path3) ### DEF1_modified ###

#Get subset of Alvinella db
bash_subset1 = get_dico_seq_subset(path1, list_query)   ### DEF2 ###

#Get subset of Paralvinella db
bash_subset2 = get_dico_seq_subset(path2, list_match)   ### DEF2 ###


## 3 ## Grab Pairwise DNA
for pair in  list_pairwise:  # list pairwise comes from '19_RBH ...'
    name1 = pair[0]    
    seq1 = bash_subset1[name1]
    name2 = pair[1]    
    seq2 = bash_subset2[name2]

    file_OUT.write(">%s\n" %name1)
    file_OUT.write("%s\n" %seq1)
    file_OUT.write(">%s\n" %name2)
    file_OUT.write("%s\n" %seq2)
    
file_OUT.close()
