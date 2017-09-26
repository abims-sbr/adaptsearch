#!/usr/bin/env python
#argv[1] : txt file with orthogroups
#argv[2] : minimal number of species to keep per group

## This script takes an output file of OrthoFinder (Orthogroups.txt), which contains a set of orthogroups,
## and rewrite it to split each orthogroup into a single fasta file.
## Beta version

""" The used output does not give the number of species per orthogroups, making filtering more difficult.
    (it's in the summar y statistics in .csv format). Nethertheless, the format of sequences IDs in the 
    AdaptSearch pipeline input files keep the record of the species : this will allow (later) the script 
    fitler_orthofinder to do the same thing than filter_fastortho.

    We could try to write a function able to read a csv file ...
"""

import os, sys, string, glob, csv
from Bio import SeqIO # BioPython

# **********************************************************************************************************************************

## PART 1 : Make a dictionary of {IDs : sequence}

""" Written for testing with the ExampleDataset of Orthofinder. Not deleted just in case.
Make a fasta copy of .faa initial files (bioconductor SeqIO does not support .faa files) """
def cpyAndRename(file):
    name = file.split('.')
    name = "%s_oneline.fasta" %name[0]
    os.system("cp %s %s" %(file, name))

""" Make the fasta copy as one sequence per line (with biopython).
Also written for testing and not deleted, just in case we need it."""
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

""" Build a hash table with gene IDs and gene sequences from fasta made from input files """
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

""" Takes a file.txt of orthogroups as parameter and return a list of list of
    orthogroups where there are at least argv[2] loci; WARNING : sequences names within 
    the groups must the same as IDs in fasta files from Filter_Assemblies. if not, the 
    dictionnary will be false. That's is why the script "format_transdecoder_headers is for."""
def formatAndFilter(orthogroups, mini):
    orthogroups = open(orthogroups, "r")
    list_orthogroups = []
    # STEP 1 - Read file into a list
    with orthogroups:
        while (1):
            group = orthogroups.readline()
            group=group[:-1] # Removes terminal '\n'
            if not group:
                break
            else:
                list_orthogroups.append(group) 
    
    # STEP 2 - Convert into a list of sublist
    list_orthogroups_format = []

    for group in list_orthogroups:
        group = string.split(group, " ")

        group.sort()
        new_group = []
        rang=-1
        for loci in group: # loop for naive paralogs filtering
            if rang == -1:
                new_group.append(loci)
                rang +=1                
            elif loci[1:3] != new_group[rang][1:3]:
                new_group.append(loci)
                rang +=1

        if len(new_group) >= mini: # Drop too small orthogroups
            list_orthogroups_format.append(new_group)
    
    return list_orthogroups_format #list_orthogroups_no_para

""" Writes each orthogroup in a fasta file. Retrieves sequences with a hash table """
def writingOutputFiles(list_orthogroups, hashTable):
    i = 1
    for group in list_orthogroups :
        length = len(group)
        name = "orthogroup_%d_%d_loci.fasta" %(i, length)
        i += 1        
        result = open(name, "w")
        with result:
            for locus in group:                
                result.write("%s\n" %locus) # write geneID. ">%s\n" before
                result.write("%s\n" %hashTable[locus]) # write sequence        

# **********************************************************************************************************************************

## PART 3 : A short summary statistics

# Return the right dimensions for a matrix
def matrixDim(listOrthogroups):
    linesNbLoci = 0
    columnsNbSpec = 9
    for group in listOrthogroups:
        if len(group) > linesNbLoci:
            linesNbLoci = len(group)
    matDim = [linesNbLoci, columnsNbSpec]
    return matDim

# Builds a matrix using the computed dimensions
def matrixConstruction(matDim):
    matrix = []
    for i in range (0,matDim[0]):
        matrix.append([0] * matDim[1])
    return matrix

# Fill the matrix (number of Orthogroups, number of sequences & species per group)
def matrixFilling(matrix, listOrthogroups):
    for group in listOrthogroups:
        listSpecs = []
        for loci in group:
            if loci[1:3] not in listSpecs:
                listSpecs.append(loci[1:3])
        matrix[len(group)-1][len(listSpecs)-1] += 1

""" print a short summary statistics. Not used anymore """
def printCounts(matrix):
    legend = []
    i=1
    for elem in matrix[0]:
        legend.append(i)
        i+=1
    
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

def writeTable(matrix, mini, filename):
    tfile=open(filename, "w")
    summary = csv.writer(tfile, delimiter='\t', quoting=csv.QUOTE_NONE)

    summary.writerow(["Details of computed orthogroups"])
    summary.writerow([])
    summary.writerow(["sequences per orthogroup", "Species per orthogroup"])


    for i in range (0, len(matrix)) :
        nb = [mini+i]     
        count = []
        for j in range (0, len(matrix[i])) :
            count.append(matrix[i][j])    
        line = nb + count    
        summary.writerow([line])

    summ = 0
    for elem in matrix:
        summ += sum(elem)  
    summary.writerow([])
    summary.writerow(["Number of orthogroups :", summ])
   
## MAIN

def main():
    print "\n-This script works on the 'Orthogroups' file output of Orthofinder to split each orthogroup in a single fasta file."
    print "-It also gets rid of orthogroups with less sequences than the number specified by the user."
    print "-It also filters (naively) paralogous genes (only one copy is left while the other are removed from groups)\n" 

    # Build hashtable
    print "  Building hashTable IDs/sequences ...\n"
    path = glob.glob('*.fasta')    
    hashTable = hashSequences(path)   

    # Open txt file with orthogroups
    print "  Reading Orthogroups.txt ..."
    print "    (Dropping orthogroups of less than " + sys.argv[2] +" loci.)\n"    
    list_orthogroups = formatAndFilter(sys.argv[1], int(sys.argv[2])) # DEF4

    # Print summary
    print "  Writing summary ...\n"
    dim = matrixDim(list_orthogroups)
    mat = matrixConstruction(dim)    
    matrixFilling(mat, list_orthogroups)    
    writeTable(mat, int(sys.argv[2]), "summary_orthogroups.csv")
    
    # Create orthogroups files
    print "  Writing output files ...\n"
    writingOutputFiles(list_orthogroups, hashTable)  

    # Move output files in a new directory
    os.system("mkdir filtered_orthogroups")
    path = glob.glob("orthogroup*")
    for file in path:
        os.system("mv %s filtered_orthogroups" %file)

    print "\t**** Results ****\n"
    printCounts(mat)
    print("\n")
    print "\t Orthogroups files are written in the directory 'filtered_orthogroups'"
    print "\t Countings are written in the file 'summary_orthogroups.csv'\n"    

if __name__ == "__main__":
    main()
