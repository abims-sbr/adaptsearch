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

####################
###### DEF 10 ######
####################
def remove_position_with_too_much_missing_data(bash_aa, bash_nuc, MIN_SPECIES_NB):

    ## 1 ## Get alignment length
    fasta_name0 = bash_aa.keys()[0]
    ln_aa = len(bash_aa[fasta_name0])

    ln_nuc = len(bash_nuc[fasta_name0])


    ## 2 ## Get positions keeped in aa alignment
    LIST_POSITION_KEEPED_aa = []
    i=0
    while i < ln_aa:
        site = []
        for fasta_name in bash_aa.keys():
            pos = bash_aa[fasta_name][i]
            #print pos

            if pos != "-" and pos != "?" and pos != "X":
                site.append(pos)
        #print site
        if len(site) >= MIN_SPECIES_NB:
            LIST_POSITION_KEEPED_aa.append(i)
        i = i+1

    ## 3 ## Get positions keeped in nuc alignment
    LIST_POSITION_KEEPED_nuc = []
    for position in LIST_POSITION_KEEPED_aa:
        position1 = position*3
        position2 = position*3 + 1
        position3 = position*3 + 2
        LIST_POSITION_KEEPED_nuc.append(position1)
        LIST_POSITION_KEEPED_nuc.append(position2)
        LIST_POSITION_KEEPED_nuc.append(position3)
        
    ## 4 ## Create entries for "filtered_bash" for aa & nuc
    filtered_bash_aa = {}
    filtered_bash_nuc = {}
    for fasta_name in bash_aa.keys():
        filtered_bash_aa[fasta_name] = ""
    for fasta_name in bash_nuc.keys():
        filtered_bash_nuc[fasta_name] = ""
    
    ## 5 ## Write "filtered_bash" for aa
    j=0
    while j < ln_aa:
        for fasta_name in bash_aa.keys():
            seq=filtered_bash_aa[fasta_name]
            pos=bash_aa[fasta_name][j]

            if j in LIST_POSITION_KEEPED_aa:
                seq = seq + pos
                filtered_bash_aa[fasta_name] = seq
        j = j + 1

    ## 6 ## Remove empty sequence
    for name in filtered_bash_aa.keys():
        seq = filtered_bash_aa[name]
        if seq == '':
            del filtered_bash_aa[name]
            
    
    ## 7 ## Write "filtered_bash" for nuc
    j=0
    while j < ln_nuc:
        for fasta_name in bash_nuc.keys():
            seq=filtered_bash_nuc[fasta_name]
            #print seq
            pos=bash_nuc[fasta_name][j]

            if j in LIST_POSITION_KEEPED_nuc:
                seq = seq + pos
                filtered_bash_nuc[fasta_name] = seq
        j = j + 1
        
    ## 8 ## Remove empty sequence
    for name in filtered_bash_nuc.keys():
        seq = filtered_bash_nuc[name]
        if seq == '':
            del filtered_bash_nuc[name]
    
    return(filtered_bash_aa, filtered_bash_nuc)
####################################


#######################
##### RUN RUN RUN #####
#######################
import string, os, time, re

### 0 ### PARAMETERS
MIN_SPECIES_NB = 4   

### 1 ### IN
path_IN1 = "./05_CDS_aa/"
#path_IN1 = "./TEST_aa"
L_IN1 = os.listdir(path_IN1)
path_IN2 = "./05_CDS_nuc/"
#path_IN2 = "./TEST_nuc"
L_IN2 = os.listdir(path_IN2)

## 2 ## OUT
path_OUT1 = "08_CDS_aa_MINIMUM_MISSING_SEQUENCES"
L_OUT1 = os.listdir(path_OUT1)
path_OUT2 = "08_CDS_nuc_MINIMUM_MISSING_SEQUENCES"
L_OUT2 = os.listdir(path_OUT2)

LOG = open("09_summary.log", "w")
LOG2 = open("08_empty_Loci_not_recorded.log", "w")

## 3 ## Clean Path OUT
if L_OUT1 != []:
    os.system("rm -fr %s/*" %path_OUT1)
if L_OUT2 != []:
    os.system("rm -fr %s/*" %path_OUT2)

## 4 ## Process files
n0=0
n1=0
n2=0
n3=0
n4=0
n5=0
n6=0
for file in L_IN1:
    print file
    file_INaa = open("%s/%s" %(path_IN1, file), "r")
    file_INnuc = open("%s/%s" %(path_IN2, file), "r")

    dico_aa = dico(file_INaa)   ### DEF 0 ###
    dico_nuc = dico(file_INnuc)   ### DEF 0 ###


    ## 4.1 ## REMOVE POSITIONS WITH TOO MUCH MISSING DATA (i.e. not enough taxa represented at each position in the alignment)
    filtered_bash_aa, filtered_bash_nuc = remove_position_with_too_much_missing_data(dico_aa, dico_nuc, MIN_SPECIES_NB)   ### DEF 10 ###
    
    ## 4.2 ## Close INPUT
    file_INaa.close()
    file_INnuc.close()

    ## 4.3 ## Change file name for output, depending the number of species remaining in the alignment
    LS = string.split(file, "_")
    ln_aa = len(filtered_bash_aa.keys())
    nb = "%dsp" %ln_aa
    new_name = LS[0] + "_" + nb + "_" + LS[2]  

    ## 4.4 ## Count and record number of sequences
    sp_nb = ln_aa

    if sp_nb == 1:
        n0=n0+1
        n1 = n1+1
    elif sp_nb == 2:
        n0=n0+1
        n2 = n2+1
    elif sp_nb == 3:
        n0=n0+1
        n3 = n3+1
    elif sp_nb == 4:
        n0=n0+1
        n4 = n4+1
    elif sp_nb == 5:
        n0=n0+1
        n5 = n5+1
    elif sp_nb == 6:
        n0=n0+1
        n6 = n6+1


    

    ## 4.5 ## Write filtered alignment in OUTPUTs
    ## aa
    if filtered_bash_aa != {}:
        OUTaa=open("%s/%s" %(path_OUT1, new_name), "w")
        for fasta_name in filtered_bash_aa.keys():
            seq_aa = filtered_bash_aa[fasta_name]
            OUTaa.write("%s\n" %fasta_name)
            OUTaa.write("%s\n" %seq_aa)
        OUTaa.close()
    else:
        LOG2.write("AA removed%s\n" %file)
        print "Alignment AA empty after filtering!!"
    # nuc
    if filtered_bash_nuc != {}:
        OUTnuc=open("%s/%s" %(path_OUT2, new_name), "w")
        for fasta_name in filtered_bash_nuc.keys():
            seq_nuc = filtered_bash_nuc[fasta_name]
            OUTnuc.write("%s\n" %fasta_name)
            OUTnuc.write("%s\n" %seq_nuc)
        OUTnuc.close() 
    else:
        LOG2.write("NUC removed%s\n\n" %file)
        print "Alignment nuc empty after filtering!!"
   
## 5 ## Print LOG
LOG.write("Total number of locus recorded  = %d\n\n" %n0)
LOG.write("Number of locus with 1 species  = %d\n" %n1)
LOG.write("Number of locus with 2 species  = %d\n" %n2)
LOG.write("Number of locus with 3 species  = %d\n" %n3)
LOG.write("Number of locus with 4 species  = %d\n" %n4)
LOG.write("Number of locus with 5 species  = %d\n" %n5)
LOG.write("Number of locus with 6 species  = %d\n" %n6)

LOG.close()
