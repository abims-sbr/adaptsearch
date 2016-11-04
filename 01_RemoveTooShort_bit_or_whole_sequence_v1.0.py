#!/usr/bin/python
## Author: Eric Fontanillas
## Last modification: 17/06/2011
## Subject: find and remove indels


###############################
##### DEF 0 : Dico fasta  #####
###############################
def dico(F2):
    #F2 = open(fasta_file_path, "r")
    dicoco = {}
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
            dicoco[fasta_name_query]=fasta_seq_query
    #F2.close()
    return(dicoco)
###################################################################################


#####################
###### DEF 1 ########
#####################
#################################
###  Create bash for genetic code
###       KEY = codon
###       VALUE = Amino Acid
#################################
def code_universel(F1):
    bash_codeUniversel = {}
    while 1:
        next = F1.readline()
        if not next: break
        L1 = string.split(next, " ")
        length1 = len(L1)
        if length1 == 3:
            key = L1[0]
            value = L1[2][:-1]
            bash_codeUniversel[key] = value
        else:
            key =  L1[0]
            value = L1[2]
            bash_codeUniversel[key] = value
    #print bash_codeUniversel
    F1.close()
    return(bash_codeUniversel)
###########################################################



#################
##### DEF 2 #####
#################
############################################################
### Test if the sequence is a multiple of 3, and if not correct the sequence to become a multiple of 3 ###
### !!!!!!!!!!!!!!!!!!!!! WEAKNESS OF THAT APPROACH = I remove extra base(s) at the end of the sequence ==> I can lost a codon, when I test ORF (as I will decay the ORF)
############################################################
def multiple3(seq):
    leng = len(seq)
    #print "\nINITIAL LENGTH = %d" %leng
    modulo = leng%3
    if modulo == 0:   # the results of dividing leng per 3 is an integer
        new_seq = seq        
    elif modulo == 1:    # means 1 extra nc (nucleotid) needs to be removed (the remaining of modulo indicate the part which is non-dividable per 3)
        new_seq = seq[:-1]  # remove the last nc
    elif modulo == 2:  # means 2 extra nc (nucleotid) needs to be removed (the remaining of modulo indicate the part which is non-dividable per 3)
        new_seq = seq[:-2]  # remove the 2 last nc
    len1 = len(new_seq)
    #print "NEW LENGTH = %d\n" %len1    
    return(new_seq, modulo)
##########################################################

###################
###### DEF 6 ######
###################
## Detect all indices corresponding to all occurance of a substring in a string
def allindices(string, sub):
    listindex=[]
    offset=0
    i = string.find(sub, offset)
    while i >= 0:
        listindex.append(i)
        i = string.find(sub, i + 1)
    return listindex
######################################################

###################
###### DEF 7 ######
###################
## Detect if methionin in the aa sequence
def detect_Methionine(seq_aa, Ortho):

    ln = len(seq_aa)

    CUTOFF_Last_50aa = ln -50

    #Ortho = 0  ## means orthologs not found
    
    ## Find all indices of occurances of "M" in a string of aa
    list_indices = allindices(seq_aa, "M")  ### DEF6 ###
         
    ## If some "M" are present, find whether the first "M" found is not in the 50 last aa (indice < CUTOFF_Last_50aa) ==> in this case: maybenot a CDS
    if list_indices != []:
        first_M = list_indices[0]
        #print first_M
        #print CUTOFF_Last_50aa
        
        if first_M < CUTOFF_Last_50aa:
            Ortho = 1  ## means orthologs found

    return(Ortho)
####################################

###################
###### DEF 8 ######
###################
## Find internal big indel
def find_internal_big_indel(seq, MIN_LENGTH_INTRON):
    ## 1 ## Built the list of sublist of consecutive gap position
    LIST = []
    sublist=[]
    ln = len(seq)
    i=0
    while i < ln:
        #if seq[i] == "-":
        if seq[i] == "?":
            sublist.append(i)  ## save gaps in sublist until a aa is found => else:
        elif seq[i] == "X":
            sublist.append(i)  ## save gaps in sublist until a aa is found => else:
        else: 
            LIST.append(sublist)  ## save the list of gap
            sublist = []  ## create new list of gap
        i = i+1
    ## if gap at the end: add the last "sublist of gap" (not done in previous loop, at it add sublist (of gaps) only when in find aa, but if gap at the end, no aa after are present, so cannot add this last sublist to the LISt of gaps
    if sublist != []:
        LIST.append(sublist)
    #print LIST

    LIST2 = []
    for element in LIST:
        if element != []:
            LIST2.append(element)

    ## 2 ## First filter of this sublist, remove sublist if stating from the bginning of the sequence or ending at the end of the sequence (i.e. respectively 5' missing data and 3' missing data)
    indice_debut = 0
    indice_fin = ln-1
    #print indice_debut
    #print indice_fin
    LIST3 = []
    for element in LIST2:
        if indice_debut in element or indice_fin in element: 
            #print element
            g = "do nothing"
        else:
            LIST3.append(element)   ## Only internal gap are recorded in LIST3 (no flanking missing data)
    #print LIST3  

    ## 3 ## Second filter of this sublist = keep as potential INTRON, the sublist contening more than MIN_LENGTH_INTRON gaps (lower ==> could be indel polymorphism)
    LIST4 = []

    for element in LIST3:
        leng = len(element)
    
        if leng > MIN_LENGTH_INTRON:
            LIST4.append(element)
    #print LIST4

    return LIST4
####################################


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
            #if seq[i] == "?":
            sublist.append(i)  ## save gaps in sublist until a aa is found => else:
        else: 
            LIST.append(sublist)  ## save the list of gap
            sublist = []  ## create new list of gap
        i = i+1
            ## if gap at the end: add the last "sublist of gap" (not done in previous loop, at it add sublist (of gaps) only when in find aa, but if gap at the end, no aa after are present, so cannot add this last sublist to the LISt of gaps
    if sublist != []:
        LIST.append(sublist)
    #print LIST

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

### 0 ### PARAMETERS
MIN_LENGTH_ALL_aa = 30 
MIN_LENGTH_BIT_OF_SEQUENCE_aa = 15
MAX_LENGTH_SMALL_INDEL = 2      ## in aa
MAX_LENGTH_SMALL_INDEL_nuc = 6  ## in nuc
MIN_SPECIES_NB = 3   

### 1 ### IN
#path_IN1 = "../../06_check_ORF_min50aa/06_CDS_with_M_aa_TEST/"
path_IN1 = "../../06_check_ORF_min50aa/06_CDS_with_M_aa/"
#path_IN1 = "./00_test_aa"
L_IN1 = os.listdir(path_IN1)
#path_IN2 = "../../06_check_ORF_min50aa/06_CDS_with_M_nuc_TEST/"
path_IN2 = "../../06_check_ORF_min50aa/06_CDS_with_M_nuc/"
#path_IN2 = "./00_test_nuc"
L_IN2 = os.listdir(path_IN2)

## 2 ## OUT
path_OUT1 = "05_CDS_aa"
L_OUT1 = os.listdir(path_OUT1)
path_OUT2 = "05_CDS_nuc"
L_OUT2 = os.listdir(path_OUT2)

LOG = open("06_summary.log", "w")

## 3 ## Clean Path OUT
if L_OUT1 != []:
    os.system("rm -fr %s/*" %path_OUT1)
if L_OUT2 != []:
    os.system("rm -fr %s/*" %path_OUT2)

## 4 ## Process files
e=0

n0 = 0
n1 = 0
n2 = 0
n3 = 0
n4 = 0
n5 = 0
n6 = 0
K=0
for file in L_IN1:
    print "PROCESS %s" %file
    file_INaa = open("%s/%s" %(path_IN1, file), "r")
    file_INnuc = open("%s/%s" %(path_IN2, file), "r")

    dico_aa = dico(file_INaa)   ### DEF 0 ###
    dico_nuc = dico(file_INnuc)   ### DEF 0 ###

    #######################
    #if K==1:
    #    print dico_nuc
    #    print dico_aa
    #    sys.exit()
    #######################
    
    
    new_bash_aa = {}
    new_bash_nuc = {}
    for fasta_name in dico_aa.keys():
        #print "\tFASTA NAME = %s" %fasta_name
        seq = dico_aa[fasta_name]
        seq_nuc = dico_nuc[fasta_name]

        if "?" in seq:
            seq = string.replace(seq, "?", "-")
        if "?" in seq_nuc:
            seq_nuc = string.replace(seq_nuc, "?", "-")
            
        #print seq
        #print len(seq)
        #print seq_nuc
        #print len(seq_nuc)

        
        ## 4.1 ## [FILTER 1] : Detect and Replace short internal indel symbole (= "-" as for other longer gaps) by a "?"

        ## aa
        list_sublist_pos = detect_short_indel(seq, MAX_LENGTH_SMALL_INDEL)   ### DEF 9 ###
        for pos_short_indels in list_sublist_pos:
            for pos in pos_short_indels:
                #print pos
                seq = seq[:pos] + "?" + seq[pos+1:]
                #print seq
        ## nuc
        list_sublist_pos = detect_short_indel(seq_nuc, MAX_LENGTH_SMALL_INDEL_nuc)   ### DEF 9 ###
        for pos_short_indels in list_sublist_pos:
            for pos in pos_short_indels:
                #print pos
                seq_nuc = seq_nuc[:pos] + "?" + seq_nuc[pos+1:]
                #print seq

        
        ## 4.2 ## [FILTER 2] : Remove short bits of sequence (<"MIN_LENGTH_BIT_OF_SEQUENCE_aa")
        LIST_sublist_aa=[]
        S1 = string.split(seq, "-")
        for element in S1:
            if len(element) > MIN_LENGTH_BIT_OF_SEQUENCE_aa:
                LIST_sublist_aa.append(element)
        #print LIST_sublist_aa

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

        
        #print new_bash_nuc

        ## 4.6 ## [FILTER 5] Remove position (site), if not at least "MIN_SPECIES_NB" species represented on this site
        ## DO THAT IN A FOLLOWING SCRIPT INTO THE PIPELINE

        ## 4.7 ## Correct the nb of sequence in the output name, if necessary
        sp_nb = len(new_bash_aa.keys())


    
        if sp_nb == 1:
            n0 = n0+1
            n1 = n1+1
        elif sp_nb == 2:
            n0 = n0+1
            n2 = n2+1
        elif sp_nb == 3:
            n0 = n0+1
            n3 = n3+1
        elif sp_nb == 4:
            n0 = n0+1
            n4 = n4+1
        elif sp_nb == 5:
            n0 = n0+1
            n5 = n5+1
        elif sp_nb == 6:
            n0 = n0+1
            n6 = n6+1
    
        lis = string.split(file, "_")

        old_nb = lis[1]
        
        nb = "%dsp" %  sp_nb

        if old_nb == nb:
            new_file = lis[0] + "_" + nb + "_" + lis[2]
        else:
            new_file = lis[0] + "_" + nb + "_NEW_" + lis[2]
            
        
    ## [FILTER 6]: print output only if at least "MIN_SPECIES_NB" species remaining in the alignment
    if sp_nb >= MIN_SPECIES_NB:
        ## 4.8 ## Open files OUT
        file_OUTaa = open("%s/%s" %(path_OUT1, new_file), "w")
        file_OUTnuc = open("%s/%s" %(path_OUT2, new_file), "w")
        ## 4.9 ## Print Output
        for fasta_name in new_bash_aa.keys():
            seq = new_bash_aa[fasta_name]
            file_OUTaa.write("%s\n" %fasta_name)
            file_OUTaa.write("%s\n" %seq)
        for fasta_name in new_bash_nuc.keys():
            seq_nuc = new_bash_nuc[fasta_name]
            file_OUTnuc.write("%s\n" %fasta_name)
            file_OUTnuc.write("%s\n" %seq_nuc)
        
        ## 4.10 ## Close Output
        file_OUTaa.close()
        file_OUTnuc.close()

    else:
        e=e+1   ## Number of locus excluded

#         ## [FILTER 6]: print output only if at least "MIN_SPECIES_NB" species remaining in the alignment
#         if sp_nb >= MIN_SPECIES_NB:
#             ## 4.8 ## Open files OUT
#             file_OUTaa = open("%s/%s" %(path_OUT1, new_file), "w")
#             file_OUTnuc = open("%s/%s" %(path_OUT2, new_file), "w")
#             ## 4.9 ## Print Output
#             for fasta_name in new_bash_aa.keys():
#                 seq = new_bash_aa[fasta_name]
#                 file_OUTaa.write("%s\n" %fasta_name)
#                 file_OUTaa.write("%s\n" %seq)
#             for fasta_name in new_bash_nuc.keys():
#                 seq_nuc = new_bash_nuc[fasta_name]
#                 file_OUTnuc.write("%s\n" %fasta_name)
#                 file_OUTnuc.write("%s\n" %seq_nuc)

#             ## 4.10 ## Close Output
#             file_OUTaa.close()
#             file_OUTnuc.close()

#         else:
#             e=e+1   ## Number of locus excluded


    ## 4.11 ## Close INPUT
    file_INaa.close()
    file_INnuc.close()

## 5 ## Print LOG
LOG.write("Total number of locus recorded  = %d\n\n" %n0)

LOG.write("Number of locus with 1 species  = %d\n" %n1)
LOG.write("Number of locus with 2 species  = %d\n" %n2)
LOG.write("Number of locus with 3 species  = %d\n" %n3)
LOG.write("Number of locus with 4 species  = %d\n" %n4)
LOG.write("Number of locus with 5 species  = %d\n\n" %n5)
LOG.write("Number of locus with 6 species  = %d\n\n" %n6)

LOG.write("Number of locus excluded (exclude if not at least %d species in the alignment)= %d\n" %(MIN_SPECIES_NB,e))

LOG.close()
