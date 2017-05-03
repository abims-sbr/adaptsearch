#!/usr/bin/env python
#coding: utf8
## AUTHOR: Eric Fontanillas
## UPDATED VERSION: 14/08/14 by Julie BAFFARD
## LAST VERSION : /01/2017 by Victor Mataigne

#MINIMUM_LENGTH = 100

# Command line : python get_locus_orthologs_part2_v2.py <zipfile:_output_pairwise>

############################
##### DEF1 : Get Pairs #####
############################
""" Read an output file from pairwise and return a list of couple [name1, sequence1, name2, sequence2]
    The purpose is to store the information of the file in a 2D tab."""
def getPairs(fasta_file_path):
    F2 = open(fasta_file_path, "r")
    list_pairwises = []
    while 1:
        next2 = F2.readline()
        if not next2:
            break
        if next2[0] == ">":
            fasta_name_query = next2[:-1] # Remove \n
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

            pairwise = [fasta_name_query,fasta_seq_query,fasta_name_match,fasta_seq_match] # List with 2 names and 2 sequences
            list_pairwises.append(pairwise) # list (2D tab : [][]) of all pairwises

    F2.close()    
    return(list_pairwises)

#################################
### DEF2 : Get a list of loci ###
#################################
""" Takes as input the list obtained in DEF1 : put in the same list all homologous sequences 
    -> If there is a couple sp1/sp2 and a couple sp1/sp3, the output is a list [sp1, sp2, sp1, sp3] """
def getListLocus(listDir):
    aListOfLocus=[]    
    
    for subfolder in listDir:
        if subfolder.startswith("19_"):
            fileIN = subfolder
            listPairwises = getPairs(fileIN) # DEF1 
            
            j = 0
            for pair in listPairwises:
                name1 = pair[0]
                name2 = pair[2]
                m = 0
            
                for locus in aListOfLocus:
                    if name1 in locus or name2 in locus:
                        locus.append(name1)
                        locus.append(name2)
                        m = 1
                if m==0:
                    newLocus = [name1, name2]
                    aListOfLocus.append(newLocus)
                
                j += 1    
    return aListOfLocus

##################################################
### DEF3 : Remove redondancy in a list of loci ###
##################################################
""" Remove the redondancy in the lists obtained with DEF2. """
def remRedondancy(aListOfLocus):

    listShortLocus = []

    for locus in aListOfLocus:
        index = aListOfLocus.index(locus)
        
        shortLocus = []
        
        for seq in locus:
            if seq not in shortLocus:
                shortLocus.append(seq)
        
        shortLocus.sort()
        
        if shortLocus not in listShortLocus:
            listShortLocus.append(shortLocus)
            
    return listShortLocus
    
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
L1 = os.listdir(PATH) # Liste des fichiers dans PATH

list_Locus = getListLocus(L1) # DEF2

list_short_locus = remRedondancy(list_Locus) # DEF3

## Control list locus length
la = len(list_Locus)
ln = len(list_short_locus)
print "\nNumber of locus before removeRedondancy = %d\n" %la
print "\nNumber of locus = %d\n" %ln

## Backup list locus with PICKLE
backup_list_LOCUS = open("02_backup_list_LOCUS", "w")
# Ecrit au format pickle dans le fichier
pickle.dump(list_short_locus, backup_list_LOCUS)
""" pickle save an object in a file and allows to load it afterwards without parsing. """

backup_list_LOCUS.close()

os.system("rm -f 19*")

