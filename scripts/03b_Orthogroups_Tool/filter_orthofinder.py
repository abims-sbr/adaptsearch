#!/usr/bin/env python
# Commandline : ./filter_orthofinder.py <Orthogroups.txt> <Nb_of_studied_species> <minimal_nb_species_per_group>

## This script takes an output file of OrthoFinder (Orthogroups.txt), which contains a set of orthogroups,
## and rewrite it to split each orthogroup into a single fasta file.
## Beta version

import os, string, glob, argparse, csv
import numpy as np
import pandas as pd

## PART 1 : Make a dictionary of {IDs : sequence}   

""" Build a hash table with gene IDs and gene sequences from fasta made from input files """
def hashSequences(path):
    hashTable = {}
    # Sequences are expected to be on one line    
    for file in path:        
        originFile = open(file, "r")
        gene = ""
        sequence = ""
        with originFile:
            while (1): # Ugly code       
                gene = originFile.readline()
                if not gene:
                    break
                gene = gene[:-1]            
                sequence = originFile.readline()
                sequence = sequence[:-1]           
                hashTable[gene] = sequence      
    return hashTable

## PART 2 : Create orthogroups file (one file per orthogroup)

""" Takes a file.txt of orthogroups as parameter and return a list of list of
    orthogroups where there are at least argv[2] loci; WARNING : sequences names within 
    the groups must the same as IDs in fasta files from Filter_Assemblies. if not, the 
    dictionnary will be false. That's is why the script "format_transdecoder_headers is for."""
def formatAndFilter(orthogroups, mini, nbspecs, hashTable, args):

    """ Builds a 2D array for a summary """
    def countings(listOrthogroups, nbcols):
        #listOrthogroups.sort().reverse()
        #nblines = len(listOrthogroups[0])
        nblines = 0    
        for group in listOrthogroups:
            if len(group) > nblines:
                nblines = len(group)
        matrix = np.array([[0]*nbcols]*nblines)
        # empty lines are avoided : first line of the frame is the line for minimal number of sequences in a group (>=mini)
        # for now, this feature diseappear when using numpy arrays and pandas :/

        for group in listOrthogroups:
            listSpecs = []
            for loci in group:
                if loci[1:3] not in listSpecs:
                    listSpecs.append(loci[1:3])
            matrix[len(group)-1][len(listSpecs)-1] += 1

        return matrix

    """ numpy 2D array in a nice dataframe """
    def asFrame(matrix) :
        index = [0]*len(matrix)
        colnames = [0]*len(matrix[0])
        index = [str(i+1)+" seqs" for i in range(len(matrix))]
        colnames = [str(i+1)+" sps" for i in range(len(matrix[0]))]
        df = pd.DataFrame(matrix, index=index, columns=colnames)
        return df # Mettre une selection pour ne renvoyer que les lignes et les colonnes qui somment > 0
        #return df.loc['4 seqs':'9 seqs'].loc[:,colnames[3:]]
        
    """ Writes each orthogroup in a fasta file. Retrieves sequences with a hash table """
    def writeOutputFile(orthogroup, hashTable, i):
        length = len(orthogroup)
        name = "orthogroup_{}_{}_loci.fasta".format(i, length)                   
        result = open(name, "w")
        with result:
            for locus in group:                
                result.write("{}\n".format(locus)) # write geneID. ">%s\n" before
                result.write("{}\n".format(hashTable[locus])) # write sequence

    # FUNCTION

    orthogroups = open(orthogroups, "r")
    list_orthogroups = []
    # STEP 1 - Read file into a list ----------------------------------------------
    with orthogroups:
        while (1):
            group = orthogroups.readline()
            group = group[:-1] # Removes terminal '\n'
            if not group:
                break
            else:
                list_orthogroups.append(group)    
    
    # STEP 2 - Paralogous filtering -----------------------------------------------
    if args.verbose :
        list_orthogroups_withpara = []
    list_orthogroups_format = []

    i = 1
    for group in list_orthogroups:
        group = string.split(group, " ") # list of lists
        group.sort()
        if args.verbose:
            if len(group) >= mini:
                list_orthogroups_withpara.append(group)            
        new_group = []
        rang=-1
        # Keep only one paralogs per species (1st encounter)
        for loci in group:
            if rang == -1:
                new_group.append(loci)
                rang +=1
            elif loci[1:3] != new_group[rang][1:3]:
                new_group.append(loci)
                rang +=1

        if len(new_group) >= mini: # Drop too small orthogroups
            list_orthogroups_format.append(new_group)
            writeOutputFile(new_group, hashTable, i)
            i += 1
    
    # STEP 3 - Print summaries ----------------------------------------------------
    if args.verbose:
        print "  Summary before paralogous filtering : \n"
        frame1 = asFrame(countings(list_orthogroups_withpara, nbspecs))
        print frame1
        #print "  Summary before paralogous filtering : \n",countings(list_orthogroups_withpara, nbspecs),"\n"
    print "  Summary after paralogous filtering : \n"
    frame2= asFrame(countings(list_orthogroups_format, nbspecs))    
    print frame2

    return list_orthogroups_format #list_orthogroups_no_para

## MAIN

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("files", help="Orthogroups.txt file from OrthoFinder")
    parser.add_argument("nbspec", type=int, help="Number of studied species")
    parser.add_argument("minspec", type=int, help="Minimal number of species to keep per group")
    parser.add_argument("-v", "--verbose", action="store_true", help="Add another summary table : countings before paralogous genes filtering")
    args = parser.parse_args()

    print "\n-This script works on the 'Orthogroups' file output of Orthofinder to split each orthogroup in a single fasta file."
    print "-It also gets rid of orthogroups with less sequences than the number specified by the user." 

    # Build hashtable
    print "  Building hashTable IDs/sequences ...\n"
    path = glob.glob('*.fasta')    
    hashTable = hashSequences(path)

    # Open txt file with orthogroups
    print "  Reading Orthogroups.txt and wrting orthogroups to separated files..."
    print "    (Dropping orthogroups of less than {} loci.)\n".format(args.minspec)
    list_orthogroups = formatAndFilter(args.files, args.minspec, args.nbspec, hashTable, args) # DEF4

    # Move output files in a new directory
    os.system("mkdir filtered_orthogroups")
    path = glob.glob("orthogroup*")
    for file in path:
        os.system("mv {} filtered_orthogroups".format(file))
    
    print "  \nOrthogroups files are written in the directory 'filtered_orthogroups'"
    print "  Countings are written in the file 'summary_orthogroups.csv'\n"    

if __name__ == "__main__":
    main()
