#!/usr/local/public/python-2.6.5/bin/python
#/usr/bin/python ==> PYTHON 2.4

## AUTHOR: Eric Fontanillas
## LAST VERSION: 14/08/14 by Julie BAFFARD


#####################################################################
##### DEF1_modified : Get Pairs(only fasta names, not squences) #####
#####################################################################
def get_pairs_modified(fasta_file_path):

    list_pairwises = []
    list_query = []
    list_match= []
    h = 0
    with open(fasta_file_path, "r") as F2:
        for name, query, name2, query2 in itertools.izip_longest(*[F2]*4):            
            if name[0] == ">":
                ## Read 1 pairwise
                fasta_name_query = name[1:-1]
                fasta_seq_query = query[:-1]
                fasta_name_match = name2[1:-1]
                fasta_seq_match = query2[:-1]
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

    print "Number of pairwises parsed = %d \n\n\n" %h
    return(list_pairwises, list_query, list_match)
##############################################


##################
###### DEF2 ######
##################
def get_dico_seq_subset(path, list_fastaNames):
    bash_subset = {}
    with open(path, "r") as F1:
        for header, seq in itertools.izip_longest(*[F1]*2):
            if header[0] == ">":
                if header[1:-1] in list_fastaNames:
                    name = header[1:-1]
                    sequence = seq[:-1]
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
os.system("cp -fr %s ./outputs_prot/ReciprocalBestHits_%s.fasta" %(path3,SHORT_FILE))

file_OUT = open("./outputs_dna/DNAalignment_corresponding_to_protein_from_RBH_%s.fasta" %SHORT_FILE, "w")

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
