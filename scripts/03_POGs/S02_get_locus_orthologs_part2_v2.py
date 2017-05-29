#!/usr/bin/env python
#coding: utf8
## AUTHOR: Eric Fontanillas
## UPDATED VERSION: 14/08/14 by Julie BAFFARD
## LAST VERSION : /01/2017 by Victor Mataigne

# Command line : python get_locus_orthologs_part2_v2.py <zipfile:_output_filter_assemblies> <integer:_minimal_number_seq_per_group> <paralogs_filtering:_yes/no>

#### DEF1 ####
## param : a list of lists of homologous loci ###
## [1rst treatment INTRA LOCUS] SELECT ONLY ORTHOLOGUOUS (list of names of orthologuous sequences) :
## test if locus is orthologuous (i.e. when a same species appear with only 1 sequence in the locus)
############################################################################################################################
def intraLocusTreatment(list_LOCUS) :
    LOCUS_without_DUPLI = []
    l = []

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

    return LOCUS_without_DUPLI

# If paralogs, keep one instead of removing entire group (update 03/2017)
def removeParalogous(list_orthogroups):
    list_orthogroups_filtered = []
    list_orthogroups.sort()

    for group in list_orthogroups:
        group.sort()
        new_group = []
        rang = -1
        for loci in group:
            if rang == -1 : # Juste pour le tout début
                new_group.append(loci)
                rang += 1
            elif loci[0:3] != new_group[rang][0:3]:
                new_group.append(loci)
                rang += 1
        list_orthogroups_filtered.append(new_group)  
    
    return list_orthogroups_filtered

### DEF2 ###
## [2nd treatment INTER LOCUS]In fact there are still some duplication in "LOCUS_without_DUPLI" 
## ==> no duplication in each loci ... OK, but several loci sharing a same sequence is still possible, 
## I will exclude this case now
############################################################################################################################
def interLocusTreatment(LOCUS_without_DUPLI) :
    list_name_seq = []
    bash_name_seq = {}

    # 3.1. Get bash = key = sequence name // value = list of locus
    for locus in LOCUS_without_DUPLI:
        for name in locus:
            if name not in list_name_seq:
                list_name_seq.append(name)
                list = [locus]
                bash_name_seq[name] = list
            else:
                bash_name_seq[name] = bash_name_seq[name] + list

    # 3.2. From the previous bash, test sequence name present in more than one loci (i.e. length of the list of loci > 1), record these paralogous loci
    k = bash_name_seq.keys()
    list_locus_paralogues = []
    for sequence in k:
        list = bash_name_seq[sequence]
        if len(list) > 1:
            for locus in list:
                list_locus_paralogues.append(locus)    ### list of list   

    # 3.3. Remove the paralogous loci from the "LOCUS_without_DUPLI" list
    for locus in list_locus_paralogues:
        if locus in LOCUS_without_DUPLI:
            LOCUS_without_DUPLI.remove(locus)

    return LOCUS_without_DUPLI

##### DEF3 : Dico fasta  #####
############################################################################################################################
# Dictionnaire {nomSequence : sequence}
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

            dicoco[fasta_name_query]=fasta_seq_query

    F2.close()
    return(dicoco)

##### DEF4 : Get sequences from list of names per locus  #####
############################################################################################################################
def get_seq_per_locus(list_locus, bash_seq, path_OUT,prefix_name):

    i=0
    for locus in list_locus:
        i = i+1
        OUT = open("%s/locus%d_%s.fasta" %(path_OUT, i, prefix_name), "w")
        
        for seq_name in locus:
            sequence = bash_seq[seq_name]
            OUT.write("%s\n" %seq_name)
            OUT.write("%s\n" %sequence)

        OUT.close()

############################################################################################################################

##### RUN RUN RUN #####
############################################################################################################################

import string, os, sys, pickle, zipfile, re

## 1 ## LOAD "PICKLE" LIST OF ALL POTENTIAL LOCUS (orthologuous and paraloguous)
############################################################################################################################
file_LOCUS = open("02_backup_list_LOCUS") # output from part1, pickled. A pickle is a Python object encoding a chain of bytes. cpickle is faster
list_LOCUS = pickle.load(file_LOCUS) # the file is depickled into a list. 
file_LOCUS.close()

L2 = [] # L2 : list of input file
zfile = zipfile.ZipFile(sys.argv[1]) # output of Filter Assemblies 
for name in zfile.namelist() :
    zfile.extract(name)
    L2.append(name)

nb=1
os.mkdir("04_LOCUS_ORTHOLOGS_UNALIGNED_perCLASS")
while nb<len(L2) :
    nb+=1
    os.mkdir("04_LOCUS_ORTHOLOGS_UNALIGNED_perCLASS/LOCUS_%i_sp" %nb)

## 2 ## intra LOCUS treatment
############################################################################################################################
LOCUS_without_DUPLI = []
if sys.argv[3] == "oui":
    LOCUS_without_DUPLI = intraLocusTreatment(list_LOCUS) #DEF1
    print "\nNUMBER OF REMAINING LOCUS AFTER INTRA LOCUS TREATMENT [REMOVE GROUPS WITH PARALOGS] = %d" %len(LOCUS_without_DUPLI)
elif sys.argv[3] == "non":
    LOCUS_without_DUPLI = removeParalogous(list_LOCUS)
    print "\nNUMBER OF REMAINING LOCUS AFTER INTRA LOCUS TREATMENT [PARALOGS NAIVE FILTERING] = %d" %len(LOCUS_without_DUPLI)

## 3 ## inter LOCUS treatment
############################################################################################################################
LOCUS_without_DUPLI2 = interLocusTreatment(LOCUS_without_DUPLI) #DEF2

print "NUMBER OF REMAINING LOCUS AFTER 2ND TREATMENT [INTER LOCUS] = %d\n\n" %len(LOCUS_without_DUPLI2)


## REMOVAL OF POGS WITH LESS THAN n SPECIES (DEFAULT IN WRAPPER : 2) ##
print "REMOVAL OF LOCUS WITH LESS THAN %s SEQUENCES\n" %sys.argv[2]

LOCUS_without_DUPLI2_FILTER_2SPECIES = []
for sublist in LOCUS_without_DUPLI2 :    
    if len(sublist) >= int(sys.argv[2]) :
        LOCUS_without_DUPLI2_FILTER_2SPECIES.append(sublist)

############################

nbr=1
dico_LOCUS_sp = {}

while nbr<len(L2) :
    nbr+=1
    dico_LOCUS_sp['%i_sp' %nbr] = []
    for locus in LOCUS_without_DUPLI2_FILTER_2SPECIES :
        if len(locus) == nbr :
            dico_LOCUS_sp['%i_sp' %nbr].append(locus)

total_number_locus = 0
list_key = []

for key in dico_LOCUS_sp.keys() :
    list_key.append(key)

list_key.sort()
list_key.reverse()
for number in list_key :    
    value = dico_LOCUS_sp[number]
    if value != [] :
        number = number.split("_")
    number = number[0]
    print "Number of species in the locus : %s" %number
    print "\tNumber of locus : %i" %len(value)
    print ""

## 4 ## GET SEQUENCES CORRESPONDING TO THE ORTHOLOGUOUS NAMES (per locus)
############################################################################################################################
BASH = {}

for subfile in L2:
    dicoco = dico(subfile)     ### DEF3 ###
    BASH.update(dicoco)        # dictionnary  {name : sequence } made from Filter assemblies output dict.update concatenate BASH and dicoco

sp = 1
while sp<len(L2) :
    sp+=1 
    get_seq_per_locus(dico_LOCUS_sp['%i_sp' %sp], BASH, "04_LOCUS_ORTHOLOGS_UNALIGNED_perCLASS/LOCUS_%i_sp" %sp, "sp%i" %sp)      ## DEF4 ##
    get_seq_per_locus(dico_LOCUS_sp['%i_sp' %sp], BASH, ".", "sp%i" %sp)        ## DEF4 ##

## 5 ## OUTPUT CONVERSION TO ZIP FORMAT
############################################################################################################################
locus_orthologs_unaligned = "^locus.*$"

f = zipfile.ZipFile("POGs_locus_orthologs_unaligned.zip", "w")

folder = os.listdir("./")

for i in folder :
    if re.match(locus_orthologs_unaligned, i) :
        f.write("./%s" %i)
