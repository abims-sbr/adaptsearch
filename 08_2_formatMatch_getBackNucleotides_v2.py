#!/usr/bin/python

## AUTHOR: Eric Fontanillas

## LAST VERSION: 26.10.2010

## DESCRIPTION: find back the nucleotide sequence corresponding to the protein sequences in file "09_onlyMatch_filtered_200bp.fasta"

## 1 ARGUMENT:  Minimum length (e.g. 100)

#MIN_LENGTH = 100

## SPLIT for string.split fasta name
SPLIT = "||"


###########
## DEF 1 ##
###########
## Generates bash, with key = fasta name; value = sequence (WITH GAP, IF ANY, REMOVED IN THIS FUNCTION)

def dico(fasta_file):
    #print "\tREDONDANT MATCHES:"
    count_fastaName=0
    F1 = open(fasta_file, "r")
    
    bash1 = {}

    w=0
    while 1:
        nextline = F1.readline()
        if not nextline:
            break        
        w = w+1



        if w%100000 == 0:
            print w
        
        if nextline[0] == ">":
            count_fastaName = count_fastaName + 1
            fasta_name = nextline[1:-1]
            nextline = F1.readline()
            sequence = nextline[:-1]
            
            #if fasta_name not in bash1.keys():
            bash1[fasta_name] = sequence
            #else:
            #    print "\t\t%s" %fasta_name
                
    print "END"
    F1.close()
    return(bash1, count_fastaName)
#####################################

###########
## DEF 2 ##
###########
## Get list of fasta names
def get_list_fasta_names(bash):
    print "Get List of fasta names"
    list_fasta_name = []

    L = bash.keys()

    print L

    for name in L:
        S1 = string.split(name, SPLIT) 
        short_name = S1[0]
        long_name = name
        L = [short_name, long_name]
        
        list_fasta_name.append(L)


    return(list_fasta_name)

###########
## DEF 3 ##
###########
## Get sequences from a list of fasta_names
def get_sequences(list_fasta_names, bash2):
    L = bash2.keys()
    ln = len(list_fasta_names)
    bash3 = {}

    i = 1
    for names in list_fasta_names:
        short_name = names[0]
        long_name = names[1]
        print "%s (%d/%d)" %(short_name, i, ln)
        
        for fasta_name in L:
            if short_name in fasta_name:
                new_fasta_name = long_name
                new_fasta_seq = bash2[fasta_name]
                bash3[new_fasta_name] = new_fasta_seq
        i = i+1
    print "END"
    return(bash3)

#####################################

#####################
#### RUN RUN RUN ####
#####################

import string, os, sys

#MIN_LENGTH = string.atoi(sys.argv[1])
WORK_DIR = sys.argv[2]

path_IN1 = "%s/09_onlyMatch_filtered.fasta" %WORK_DIR
path_IN2 = sys.argv[1]  ## Initially the DB in blast input

path_OUT = "%s/09_onlyMatch_filtered_nucleotideBACK.fasta" %WORK_DIR
file_OUT = open(path_OUT, "w")

print "\n****** BUILT DICTIONNARY (09_OnlyMatch...) ****** \n"
bash1, count_FastaName1 = dico(path_IN1)                              ### DEF1 ###
ln1 = len(bash1.keys())
print "%d/%d" %(ln1, count_FastaName1)
#cPickle.dump(bash1, open("dico1.file", 'wb'))
#bash1 = cPickle.load(open("dico1.file"))


print "\n****** BUILT DICTIONNARY (02_CFOF_...) ****** \n"
bash2, count_FastaName2  = dico(path_IN2)                              ### DEF1 ###
ln2 = len(bash2.keys())
print "%d/%d" %(ln2, count_FastaName2)
#cPickle.dump(bash2, open("dico2.file", 'wb'))
#bash2 = cPickle.load(open("dico2.file"))


list_fasta_names = get_list_fasta_names(bash1)      ### DEF2 ###
#cPickle.dump(bash2, open("list1.file", 'wb'))
#list_fasta_names = cPickle.load(open("list1.file"))


#print list_fasta_names[1:100]

ln = len(list_fasta_names)
#print ln
#print list_fasta_names


print "get corresponding sequences ..."
bash3 = get_sequences(list_fasta_names, bash2)      ### DEF3 ###
print "sequences grabbed!!"

print "Print output ..."
## Print OUTPUT
for fasta_name in bash3:    
    file_OUT.write(">%s\n" %fasta_name)
    file_OUT.write("%s\n" %bash3[fasta_name])    


file_OUT.close()
