#!/usr/bin/env python
#argv[1] : txt file with orthogroups
#argv[2] : minimal number of species to keep per group

## This script takes an output file of OrthoFinder (Orthogroups.txt), which contains a set of orthogroups,
## and rewrite it to split each orthogroup into a single fasta file.
## It does not remove paralogous yet
## Alpha version

""" The used output does not give the number of species per orthogroups, making filtering more difficult.
    (it's in the summar y statistics in .csv format). Nethertheless, the format of sequences IDs in the 
    AdaptSearch pipeline input files keep the record of the species : this will allow (later) the script 
    fitler_orthofinder to do the same thing than filter_fastortho.

    We could try to write a function able to read a csv file ...
"""

import os, sys, string, glob
from Bio import SeqIO # BioPython

# **********************************************************************************************************************************

## PART 1 : Make a dictionary of IDs : sequence

""" DEF1 : make a fasta copy of .faa initial files (bioconductor SeqIO does not support .faa files) """
def cpyAndRename(file):
    name = file.split('.')
    name = "%s_oneline.fasta" %name[0]
    os.system("cp %s %s" %(file, name))

""" DEF 2 :make the fasta copy as one sequence per line (with biopython)"""
def seqOneLine(file):
    n = file.split('.')    
    name = "%s_oneline.fasta" %n[0]   
    
    new_file = open(name, "w")
    with new_file:
        for seq_record in SeqIO.parse(file, "fasta"):
            gid = seq_record.id
            gid = gid.split("(")
            gid = gid[0]
            new_file.write(gid)
            new_file.write("\n")
            new_file.write(str(seq_record.seq))
            new_file.write("\n")    

""" DEF 3 :Build a hash table with gene IDs and gene sequences from fasta made from input files (with DEF1 and DEF2)"""
def hashSequences(path):
    hashTable = {}
    # Sequences are expected to be on one line    
    for file in path:        
        originFile = open(file, "r")
        gene = ""
        sequence = ""
        with originFile:
            while (1):           
                gene = originFile.readline()
                if not gene:
                    break
                gene = gene[:-1]            
                sequence = originFile.readline()
                sequence = sequence[:-1]           
                hashTable[gene] = sequence
        #originFile.close() #no need if with is used       
    return hashTable

# **********************************************************************************************************************************

## PART 2 : Create orthogroups file (one file per orthogroup)

""" remove paralogous from orthogroups (called during DEF6, AFTER DEFY )"""
def removeParalogous(list_orthogroups):
    list_orthogroups_filtered = []
    list_orthogroups.sort()

    for group in list_orthogroups:
        group.sort()
        new_group = []
        rang = -1
        for loci in group:
            if rang == -1 :
                new_group.append(loci)
                rang += 1
            elif loci[0:2] != new_group[rang][0:2]:
                new_group.append(loci)
                rang += 1
        list_orthogroups_filtered.append(new_group)  
    
    return list_orthogroups_filtered

""" DEF4 : takes a file.txt of orthogroups as parameter and return a list of list of
    orthogroups where there are at least argv[2] loci """
def formatAndFilter(orthogroups, min):
    orthogroups = open(orthogroups, "r")
    list_orthogroups = []
    # STEP 1 - Read file into a list
    with orthogroups:
        while (1):
            group = orthogroups.readline()
            if not group:
                break
            else:
                list_orthogroups.append(group)    

    # STEP 2 - Convert into a list of sublist
    list_orthogroups_format = []

    for group in list_orthogroups:
        group = string.split(group, " ")
        group = group[1:]
        list_orthogroups_format.append(group)    
    
    # STEP 4 - Remove orthogroups < argv[2] (less than argv|2) loci)
    list_orthogroups_filtered = []
    [list_orthogroups_filtered.append(group) for group in list_orthogroups_format if len(group) >= min]

    # STEP 5 - Keep only one paralog if any :
    list_orthogroups_no_para = removeParalogous(list_orthogroups_filtered)
    
    # Remove terminal '\n' in the last id
    """
    for group in list_orthogroups_filtered:
        if group[-1][-1] == '\n':
            group[-1] = group[-1][0:-1]
    
    return list_orthogroups_filtered
    """
    for group in list_orthogroups_no_para:
        if group[-1][-1] == '\n':
            group[-1] = group[-1][0:-1]
    
    return list_orthogroups_no_para

""" DEF 5 : writes each orthogroup in a fasta file. Retrieves sequences with a hash table """
def writingOutputFiles(list_orthogroups, hashTable):
    i = 1
    for group in list_orthogroups :
        length = len(group)
        name = "orthogroup_%d_%d_loci.fasta" %(i, length)
        i += 1        
        result = open(name, "w")
        with result:
            for locus in group:                
                result.write(">%s\n" %locus) # write geneID
                result.write("%s\n" %hashTable[locus]) # write sequence        

# **********************************************************************************************************************************

## PART 3 : A short summary statistics printed on the screen

""" DEF 6 : print a (very) short summary statistics""" 

def matrixDim(listOrthogroups):
    linesNbLoci = 0
    columnsNbSpec = 9
    for group in listOrthogroups:
        if len(group) > linesNbLoci:
            linesNbLoci = len(group)
    matDim = [linesNbLoci, columnsNbSpec]
    return matDim

def matrixConstruction(matDim):
    matrix = []
    for i in range (0,matDim[0]):
        matrix.append([0] * matDim[1])
    return matrix

def matrixFilling(matrix, listOrthogroups):
    for group in listOrthogroups:
        listSpecs = []
        for loci in group:
            if loci[0:2] not in listSpecs:
                listSpecs.append(loci[0:2])
        matrix[len(group)-1][len(listSpecs)-1] += 1

def printCounts(matrix, listOrthogroups):
    legend = []
    i=1
    for elem in matrix[0]:
        legend.append(i)
        i+=1

    print "\nThis script works on the 'Orthogroups' file output of Orthofinder to split each orthogroup\nin a single fasta file."
    print "\nThe script get rids of orthogroups with less sequences than the number specified by the user."        
    
    print "\n\t**** Results ****\n"

    print "\tSpecies per orthogroup : |",
    for elem in legend:
        print elem,   
    print "\n\t-------------------------|---------------------------"

    i=1 # int(sys.argv[3])
    for a in matrix:
        if sum(a) !=0 : # N'affiche pas les lignes ou il n'y a rien
            print "\t loci per orthogroup :",i, "|",
            for elem in a:
                print elem,
            print ""         
        i +=1 

    print "\nTotal of orthogroups :", len(listOrthogroups)
    print "\nNot recorded : orthogroups of less than " + sys.argv[2] +" loci.\n"  

# **********************************************************************************************************************************  
   
## MAIN

def main():
    # Build hashtable    
    path = glob.glob('*.fasta')
    for file in path:
        #cpyAndRename(file) # DEF1
        seqOneLine(file) # DEF2
    path = glob.glob('*_oneline.fasta')
    hashTable = hashSequences(path) # DEF3    

    # Open txt file with orthogroups    
    list_orthogroups = formatAndFilter(sys.argv[1], int(sys.argv[2])) # DEF4
    
    # Print summary
    dim = matrixDim(list_orthogroups)
    mat = matrixConstruction(dim)    
    matrixFilling(mat, list_orthogroups)    
    printCounts(mat, list_orthogroups) 
    
    # Create orthgroups files
    writingOutputFiles(list_orthogroups, hashTable) # DEF5    

    # Move output files in a new directory
    os.system("mkdir filtered_orthogroups")
    path = glob.glob("orthogroup*")
    for file in path:
        os.system("mv %s filtered_orthogroups" %file)
    os.system("rm *oneline.fasta")

if __name__ == "__main__":
    main()
