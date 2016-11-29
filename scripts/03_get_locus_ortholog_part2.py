#!/usr/bin/python

##############################
##### DEF1 : Dico fasta  #####
##############################
def dico(fasta_file_path):
    F2 = open(fasta_file_path, "r")
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
#             next3 = F2.readline()    ## jump one empty line (if any after the sequence)
#             #next3 = F2.readline()
#             fasta_name_match = next3[:-1]
#             Sn = string.split(fasta_name_match, "||")
#             fasta_name_match = Sn[0]

#             next3 = F2.readline()
#             fasta_seq_match = next3[:-1]

            #pairwise = [fasta_name_query,fasta_seq_query,fasta_name_match,fasta_seq_match]
            dicoco[fasta_name_query]=fasta_seq_query
            #dicoco[fasta_name_match]=fasta_seq_match

            
            ## ADD pairwise with condition
            #list_pairwises.append(pairwise)
    F2.close()
    return(dicoco)
###################################################################################

##############################################################
##### DEF2 : Get sequences from list of names per locus  #####
##############################################################
def get_seq_per_locus(list_locus, bash_seq, path_OUT,prefix_name):

    print "\nGet seq per locus\n\n"

    print "\tLIST LOCUS = %s\n\n" %list_locus
    
    print "\tBASH SEQ KEYS = %s\n\n" %bash_seq.keys()    
    
    

    i=0
    for locus in list_locus:
        i = i+1
        path_locus_OUT = "%s/%s_%d.fasta" %(path_OUT, prefix_name, i)
        OUT = open(path_locus_OUT, "w")
        
        for seq_name in locus:
            sequence = bash_seq[seq_name]
            OUT.write("%s\n" %seq_name)
            OUT.write("%s\n" %sequence)

        OUT.close()
###################################################################################



#######################
##### RUN RUN RUN #####
#######################
import string, os, sys, pickle

## 1. LOAD "PICKLE" LIST OF ALL POTENTIAL LOCUS (orthologuous and paraloguous)
## 2. SELECT ONLY ORTHOLOGUOUS (list of names of orthologuous sequences)
## 3. GET SEQUENCES CORRESPONDING TO THE ORTHOLOGUOUS NAMES (per locus)


############################################################
## 1. LOAD "PICKLE" LIST OF ALL POTENTIAL LOCUS
############################################################
file_LOCUS = open("02_backup_list_LOCUS")
list_LOCUS = pickle.load(file_LOCUS) 
file_LOCUS.close()

LOG = open("04_ResultsParaloguesRemoval.log", "w")

######################################################
## 2. [1rst treatment INTRA LOCUS] SELECT ONLY ORTHOLOGUOUS: test if locus is orthologuous (i.e. when a same species appear with only 1 sequence in the locus)
######################################################
LOCUS_without_DUPLI = []
for locus in list_LOCUS:
    ## remove redondancy (i.e. same names present several times as it is expected)
    l = []
    for name in locus:
        if name not in l:
            if name[-2:] == "_S":
                name = name[:-2]
            l.append(name)

    ## Now test if only one name per species, if not ==> duplication
    list_initials = []
    D=0
    for name in l:
        initials = name[1:3]

        if initials not in list_initials:
            list_initials.append(initials)            
        else:  # DUPLICATION DETECTED
            D=1

    if D==0:  # means no duplication detected
        LOCUS_without_DUPLI.append(l)

LOG.write("NUMBER OF REMAINING LOCUS AFTER 1RST TREATMENT [INTRA LOCUS] = %d\n" %len(LOCUS_without_DUPLI))

print LOCUS_without_DUPLI[1:10]
################################################################
## 3. [2nd treatment INTER LOCUS] In fact there are still some duplication in "LOCUS_without_DUPLI" ==> no duplication in each loci ... OK, but several loci sharing a same sequence is still possible, I will exclude this case now
################################################################
list_name_seq = []
bash_name_seq = {}

# 3.1. Get bash = key = sequence name // value = list of locus
for locus in LOCUS_without_DUPLI:
    for sequence in locus:
        if sequence not in list_name_seq:
            list_name_seq.append(sequence)
            list = [locus]
            bash_name_seq[sequence] = list

        else:
            bash_name_seq[sequence] = bash_name_seq[sequence] + list

# 3.2. From the previous bash, test sequence name present in more than one loci (i.e. length of the list of loci > 1), record these paralogous loci
k = bash_name_seq.keys()
list_locus_paralogues = []
for sequence in k:
    list = bash_name_seq[sequence]
    if len(list) > 1:
        for locus in list:
            list_locus_paralogues.append(locus)    ### list of list

print list_locus_paralogues[1:20]   

# 3.3. Remove the paralogous loci from the "LOCUS_without_DUPLI" list

for locus in list_locus_paralogues:
    #print "\t%s" %locus
    if locus in LOCUS_without_DUPLI:
        LOCUS_without_DUPLI.remove(locus)


LOG.write("NUMBER OF REMAINING LOCUS AFTER 2ND TREATMENT [INTER LOCUS] = %d\n\n" %len(LOCUS_without_DUPLI))

##################################        

list_LOCUS_2sp = []
list_LOCUS_3sp = []
list_LOCUS_4sp = []
list_LOCUS_5sp = []
list_LOCUS_6sp = []

for locus in LOCUS_without_DUPLI:
    if len(locus) == 2:
        list_LOCUS_2sp.append(locus)
    if len(locus) == 3:
        list_LOCUS_3sp.append(locus)
    if len(locus) == 4:
        list_LOCUS_4sp.append(locus)
    if len(locus) == 5:
        list_LOCUS_5sp.append(locus)
    if len(locus) == 6:
        list_LOCUS_6sp.append(locus)

LOG.write("Number of locus with 2 species: %d\n" %len(list_LOCUS_2sp))
LOG.write("Number of locus with 3 species: %d\n" %len(list_LOCUS_3sp))
LOG.write("Number of locus with 4 species: %d\n" %len(list_LOCUS_4sp))
LOG.write("Number of locus with 5 species: %d\n\n" %len(list_LOCUS_5sp))
LOG.write("Number of locus with 6 species: %d\n\n" %len(list_LOCUS_6sp))

#########################################################################
## 4. GET SEQUENCES CORRESPONDING TO THE ORTHOLOGUOUS NAMES (per locus)
#########################################################################

path_TRANSCRIPTS = "../../tmp/01_formated_INPUT/"
L2 = os.listdir(path_TRANSCRIPTS)

BASH = {}

for subfile in L2:
    fasta_file_path = "%s/%s" %(path_TRANSCRIPTS, subfile)
    dicoco = dico(fasta_file_path)

    BASH.update(dicoco)

LOG.write("TOTAL NUMBER OF LOCI PROCESSED = %d\n" %len(BASH.keys()))
#print BASH.keys()[1:10]

## list_LOCUS_2sp
get_seq_per_locus(list_LOCUS_2sp, BASH, "04_LOCUS_ORTHOLOGS_UNALIGNED_perCLASS/LOCUS_2sp","locus_2sp")   ### DEF 2 ###
get_seq_per_locus(list_LOCUS_2sp, BASH, "04_LOCUS_ORTHOLOGS_UNALIGNED_ALL","locus_2sp") ### DEF 2 ###


## list_LOCUS_3sp
get_seq_per_locus(list_LOCUS_3sp, BASH, "04_LOCUS_ORTHOLOGS_UNALIGNED_perCLASS/LOCUS_3sp","locus_3sp")   ### DEF 2 ###
get_seq_per_locus(list_LOCUS_3sp, BASH, "04_LOCUS_ORTHOLOGS_UNALIGNED_ALL","locus_3sp")   ### DEF 2 ###


## list_LOCUS_4sp
get_seq_per_locus(list_LOCUS_4sp, BASH, "04_LOCUS_ORTHOLOGS_UNALIGNED_perCLASS/LOCUS_4sp","locus_4sp")   ### DEF 2 ###
get_seq_per_locus(list_LOCUS_4sp, BASH, "04_LOCUS_ORTHOLOGS_UNALIGNED_ALL","locus_4sp")   ### DEF 2 ###


## list_LOCUS_5sp
get_seq_per_locus(list_LOCUS_5sp, BASH, "04_LOCUS_ORTHOLOGS_UNALIGNED_perCLASS/LOCUS_5sp","locus_5sp")   ### DEF 2 ###
get_seq_per_locus(list_LOCUS_5sp, BASH, "04_LOCUS_ORTHOLOGS_UNALIGNED_ALL","locus_5sp")   ### DEF 2 ###


## list_LOCUS_6sp
get_seq_per_locus(list_LOCUS_6sp, BASH, "04_LOCUS_ORTHOLOGS_UNALIGNED_perCLASS/LOCUS_6sp","locus_6sp")   ### DEF 2 ###
get_seq_per_locus(list_LOCUS_6sp, BASH, "04_LOCUS_ORTHOLOGS_UNALIGNED_ALL","locus_6sp")   ### DEF 2 ###
