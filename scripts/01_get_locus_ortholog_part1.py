#!/usr/bin/python
## AUTHOR: Eric Fontanillas
## LAST VERSION: 14/08/14 by Julie BAFFARD

#MINIMUM_LENGTH = 100


############################
##### DEF1 : Get Pairs #####
############################
def get_pairs(fasta_file_path):
    F2 = open(fasta_file_path, "r")
    list_pairwises = []
    while 1:
        next2 = F2.readline()
        if not next2:
            break
        if next2[0] == ">":
            fasta_name_query = next2[:-1]
            Sn = string.split(fasta_name_query, "||")
            fasta_name_query = Sn[0]
            next3 = F2.readline()
            fasta_seq_query = next3[:-1]
            next3 = F2.readline()    ## jump one empty line (if any after the sequence)
            fasta_name_match = next3[:-1]
            Sn = string.split(fasta_name_match, "||")
            fasta_name_match = Sn[0]

            next3 = F2.readline()
            fasta_seq_match = next3[:-1]

            pairwise = [fasta_name_query,fasta_seq_query,fasta_name_match,fasta_seq_match]
            list_pairwises.append(pairwise)

    F2.close()
    return(list_pairwises)
##############################################


#######################
##### RUN RUN RUN #####
#######################
import string, os, time, re, sys, pickle, zipfile #modif julie : ajout de zipfile

## 1 ## INPUT/OUTPUT
#extraction du fichier zip contenant les sorties
zfile = zipfile.ZipFile(sys.argv[1])
for name in zfile.namelist() :
	zfile.extract(name)

PATH = "./"
L1 = os.listdir(PATH)

## 2 ## RUN
##Get list locus
list_LOCUS=[]
for subfolder in L1:
    if subfolder.startswith("19_") : 
    	fileIN = subfolder
	list_pairwises = get_pairs(fileIN)          ### DEF1 ###

    	jj=0
    	for pair in list_pairwises:
            name1 = pair[0]
            name2 = pair[2]         
            m=0

            ## If one of the 2 names yet present ==> complete the locus
            for locus in list_LOCUS:
            	if name1 in locus or name2 in locus:
                    locus.append(name1)
                    locus.append(name2)
                    m=1

            ## If no names was yet present in locus ==> create the locus
            if m==0:
                new_locus = [name1, name2]
            	list_LOCUS.append(new_locus)

            ## Compteur
            jj = jj+1

## Remove redondancy in list_LOCUS : Several locus in the list_locus, should be the same (but not in the same order), depending in which order they were created, we should remove them

list_short_LOCUS = []
for locus in list_LOCUS:
    index = list_LOCUS.index(locus)

    short_locus = []

    ## Remove redondant sequence name
    for seq in locus:
        if seq not in short_locus:
            short_locus.append(seq)
            
    ## Sort the list of sequence name (by alphabetical order)
    short_locus.sort()

    ## Add the locus to the new list (list_short_LOCUS) if not yet present (avoid redondant locus)
    if short_locus not in list_short_LOCUS:
        list_short_LOCUS.append(short_locus)
    
list_LOCUS = list_short_LOCUS

## Control list locus length
ln = len(list_LOCUS)
print "\nNumber of locus = %d\n" %ln

## Backup list locus with PICKLE
backup_list_LOCUS = open("02_backup_list_LOCUS", "w")
pickle.dump(list_LOCUS, backup_list_LOCUS) 

backup_list_LOCUS.close()
