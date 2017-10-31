#!/usr/bin/python
## Author: Eric Fontanillas
## Last modification: 03/09/14 by Julie BAFFARD

## Description : find and remove indels


###################
###### DEF 9 ######
###################
def detect_short_indel(seq,MAX_LENGTH_SMALL_INDEL):
    ## 1 ## Built the list of sublist of consecutive gap position
    LIST = []
    sublist=[]
    ln = len(seq)
    i=0
    while i < ln:
        if seq[i] == "-":
            sublist.append(i)  ## save gaps in sublist until a aa is found => else:
        else: 
            LIST.append(sublist)  ## save the list of gap
            sublist = []  ## create new list of gap
        i = i+1
            ## if gap at the end: add the last "sublist of gap" (not done in previous loop, at it add sublist (of gaps) only when in find aa, but if gap at the end, no aa after are present, so cannot add this last sublist to the LISt of gaps
    if sublist != []:
        LIST.append(sublist)

    ## 2 ## keep only the records of the small indel (<MAX_LENGTH_SMALL_INDEL)
    list_of_sublist_positions = []
    for element in LIST:
        if element != [] and len(element)<=MAX_LENGTH_SMALL_INDEL:
            list_of_sublist_positions.append(element)

    return(list_of_sublist_positions)
####################################




#######################
##### RUN RUN RUN #####
#######################
import string, os, time, re, sys
from dico import dico

### 0 ### PARAMETERS
MIN_LENGTH_ALL_aa = int(sys.argv[3])-20
MIN_LENGTH_BIT_OF_SEQUENCE_aa = int(sys.argv[4])
MAX_LENGTH_SMALL_INDEL = 2      ## in aa
MAX_LENGTH_SMALL_INDEL_nuc = 6  ## in nuc
MIN_SPECIES_NB = int(sys.argv[1])
MAX_sp = MIN_SPECIES_NB
dicoco = {}
dico_dico = {}
list_sp = []
list_new_file = []
n0 = 0
e=0
j=0
i=1

### 1 ### IN
if sys.argv[2] == "oui" :
    path_IN1 = "./06_CDS_with_M_aa/"
    L_IN1 = os.listdir(path_IN1)
    path_IN2 = "./06_CDS_with_M_nuc/"
    L_IN2 = os.listdir(path_IN2)
elif sys.argv[2] == "non" :
    path_IN1 = "./05_CDS_aa/"
    L_IN1 = os.listdir(path_IN1)
    path_IN2 = "./05_CDS_nuc/"
    L_IN2 = os.listdir(path_IN2)    

## 2 ## OUT
os.mkdir("07_CDS_aa")
path_OUT1 = "07_CDS_aa"
os.mkdir("07_CDS_nuc")
path_OUT2 = "07_CDS_nuc"

for file in L_IN1:
    file_INaa = open("%s/%s" %(path_IN1, file), "r")
    file_INnuc = open("%s/%s" %(path_IN2, file), "r")

    dico_aa = dico(file_INaa)   ### DEF 0 ###
    dico_nuc = dico(file_INnuc)   ### DEF 0 ###
       
    new_bash_aa = {}
    new_bash_nuc = {}
    for fasta_name in dico_aa.keys():
        seq = dico_aa[fasta_name]
        seq_nuc = dico_nuc[fasta_name]

        if "?" in seq:
            seq = string.replace(seq, "?", "-")
        if "?" in seq_nuc:
            seq_nuc = string.replace(seq_nuc, "?", "-")

        
        ## 4.1 ## [FILTER 1] : Detect and Replace short internal indel symbole (= "-" as for other longer gaps) by a "?"
        ## aa
        list_sublist_pos = detect_short_indel(seq, MAX_LENGTH_SMALL_INDEL)   ### DEF 9 ###
        for pos_short_indels in list_sublist_pos:
            for pos in pos_short_indels:
                seq = seq[:pos] + "?" + seq[pos+1:]
        ## nuc
        list_sublist_pos = detect_short_indel(seq_nuc, MAX_LENGTH_SMALL_INDEL_nuc)   ### DEF 9 ###
        for pos_short_indels in list_sublist_pos:
            for pos in pos_short_indels:
                seq_nuc = seq_nuc[:pos] + "?" + seq_nuc[pos+1:]

        
        ## 4.2 ## [FILTER 2] : Remove short bits of sequence (<"MIN_LENGTH_BIT_OF_SEQUENCE_aa")
        LIST_sublist_aa=[]
        S1 = string.split(seq, "-")
        for element in S1:
            if len(element) > MIN_LENGTH_BIT_OF_SEQUENCE_aa:
                LIST_sublist_aa.append(element)


        ## 4.3 ## [FILTER 3] : Remove all the sequence if the total length of all subsequences < "MIN_LENGTH_ALL_aa")
        seq_all = ""
        for bit_of_sequence in LIST_sublist_aa:
            seq_all = seq_all + bit_of_sequence

        if len(seq_all) < MIN_LENGTH_ALL_aa:
            LIST_sublist_aa = []


        ## 4.4 ## [FILTER 4] : Detect sublist position in the original sequence, and recreate the filtered sequence from these positions:
        seq_gap = "-" * len(seq)    ## 4.4.1 ## generate a sequence with only gaps inside
        seq_gap_nuc = "-" * len(seq_nuc)
        
        for subsequence in LIST_sublist_aa:
            ## aa
            START = string.find(seq, subsequence)
            END = START + len(subsequence)
            seq_gap = seq_gap[:START] + seq[START:END] + seq_gap[END:]  ## 4.4.2 ## and then replace the correponding gaps by coding subsequence in the sequence
            ## nuc
            START_nuc = START*3
            END_nuc = END*3
            seq_gap_nuc = seq_gap_nuc[:START_nuc] + seq_nuc[START_nuc:END_nuc] + seq_gap_nuc[END_nuc:]

        
        ## 4.5 ## Save new sequence in bash if not empty
        seq_empty_test = string.replace(seq_gap, "-", "")
        if seq_empty_test != "":
            new_bash_aa[fasta_name] = seq_gap

        seq_empty_test = string.replace(seq_gap_nuc, "-", "")
        if seq_empty_test != "":
            new_bash_nuc[fasta_name] = seq_gap_nuc

    ## 4.6 ## Correct the nb of sequence in the output name, if necessary
    sp_nb = len(new_bash_aa.keys())
    lis = string.split(file, "_")
    new_lis = string.split(lis[1], ".")       
    nb = "sp%d" %  sp_nb

    old_nb = new_lis[0]
    max_old_nb = string.split(old_nb, "sp")
    max_old_nb = "".join(max_old_nb)

    if old_nb == nb:
        new_file = lis[0] + "_" + nb + "." + new_lis[1]
    else:
        new_file = lis[0] + "_NEW_"+ nb + "." + new_lis[1]
    list_new_file.append(new_file)
    dico_dico[new_file] = [new_bash_aa, new_bash_nuc]
    n0+=1 #number of locus traited

    if list_sp == [] and nb!="sp0" :
        list_sp.append(nb)
    else :
	if nb not in list_sp and nb != "sp0" :
	    list_sp.append(nb)
 
    file_INaa.close()
    file_INnuc.close()
 
# [FILTER 5]: check if the number of locus with the max number of species isn't 0
#if it is : MIN_SPECIES_NB - 1
if len(list_sp) < MIN_SPECIES_NB :
    MIN_SPECIES_NB = len(list_sp)

## [FILTER 6]: print output only if at least "MIN_SPECIES_NB" species remaining in the alignment
for name in list_new_file :
    dicoo = dico_dico[name]
    dico_aa = dicoo[0]
    dico_nuc = dicoo[1]
    sp_nbre = len(dico_aa.keys())

    if sp_nbre >= MIN_SPECIES_NB :
        file_OUTaa = open("%s/%s" %(path_OUT1, name), "w")
        file_OUTnuc = open("%s/%s" %(path_OUT2, name), "w")

        for fasta_name in dico_aa.keys() :
	    seq_aa = dico_aa[fasta_name]
            file_OUTaa.write("%s\n" %fasta_name)
            file_OUTaa.write("%s\n" %seq_aa)
	for fasta_name in dico_nuc.keys() :
	    seq_nuc = dico_nuc[fasta_name]
            file_OUTnuc.write("%s\n" %fasta_name)
            file_OUTnuc.write("%s\n" %seq_nuc)

        file_OUTaa.close()
        file_OUTnuc.close()

    else:
        e+=1


###Print
if sys.argv[2] == "oui" :
    print "\nIn locus with CDS considering Methionine : \n"
else :
    print "\nIn locus with CDS regardless of the Methionine : \n"

print "*************** 1st filter : selection of the locus ***************"
print "\nTotal number of locus recorded  = %d" %n0 

list_sp.sort()
for sp in list_sp :
    dicoco[sp] = []
    for files in list_new_file :
        nb_sp = string.split(files, "_NEW")
	nb_sp = "".join(nb_sp)
	nb_sp = string.split(nb_sp,'_')
    	nb_sp = nb_sp[1]
        nb_sp = string.split(nb_sp, ".")
        nb_sp = nb_sp[0]
    	if nb_sp == sp :
	    dicoco[sp].append(files)
    new_sp = sp.split("sp")
    new_sp = int(new_sp[1])
    if new_sp == i :
        print "\tNumber of locus with %d species : %d" %(i, len(dicoco[sp]))
    elif i <= MAX_sp :
        print "\tNumber of locus with %d species : 0" %i 
    i+=1
if len(list_sp) != int(MAX_sp) :
    print "\tNumber of locus with %d species : 0" %int(MAX_sp)

print "Number of locus excluded (exclude if not at least %d species in the alignment)= %d\n" %(MIN_SPECIES_NB,e)
