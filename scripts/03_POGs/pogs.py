#!/usr/bin/env python
# coding: utf8
# September 2017 - Author : Victor Mataigne (Station Biologique de Roscoff - ABiMS)
# Command line : ./pogsPOO.py <list_of_input_files_separated_by_commas> <minimal number of species per group> <total_number_of_studied_species> [-v) [-p]

"""
What it does:
    - pogs.py parses output files from the "pairwise" tool of the AdaptSearch suite and proceeds to gather genes in orthogroups (using transitivity).
    - A minimal number of species per group has to be set.

BETA VERSION :

Errors corrected / spotted (need to test more to be sure) :

    - In makeOrthogroups()

        - ERROR 1

            - Some groups can be splitted :
            During the llop, if we have an orthogroup like this :
                (seq1, seq2, seq3, seq4)
            And two pairs to analyze like this :
                (seq5, seq6)
                (seq4, seq6)
            Then the the pair (seq5, seq6) won't be added to the orthogroup if it is parsed before the pair (seq4, seq6)
            
            - The problem is solved by running the loop a second time, but we need more testing to know if two loops are enought for all cases 
              or if we have to build a stronger code.
        
        - ERROR 2
            - With one reduced datatest, the 'intersection' native set method did not work when comparing pairs with two identical headers. 
            - Strangely, the same comparison did not happen with a bigger dataset
"""

import os, argparse
import numpy as np
import pandas as pd

""" Definition of a locus : header + sequence + a tag """
class Locus:    

    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence
        self.tagged = False

    def __str__(self):
        return "{}{}".format(self.header, self.sequence)

    def __eq__(self, other):
        return self.getHeader() == other.getHeader()  # Test if two loci are the same

    def __hash__(self):
        return hash((self.header, self.sequence))  # Make the object iterable and hashable

    def getHeader(self):
        return self.header

    def getSequence(self):
        return self.sequence

    def getTag(self):
        return self.tagged

    def prettyPrint(self): 
        # Used for debugging : print "{ Header : ", self.header[0:-1], "Tag : ", self.tagged, " }"
        print "[ Header : {header} ]".format(header=self.header[0:-1])

    def prettyPrint2(self):
        print "[ Header : {header} Sequence : {sequence} ]".format(header=self.header[0:-1], sequence=self.sequence[0:-1])

""" Applies the getPairwiseCouple() function to a list of files and return a big list with ALL pairwises couples 
    Returns a list of sets (2 items per set) """
def getListPairwiseAll(listPairwiseFiles):    
    
    # Sub-Function

    """ Reads an output file from the 'Pairwise' tool (AdaptSearch suite) and returns its content into a list 
        Returns a list of sets (2 items per set) """
    def getPairwiseCouple(pairwiseFile):        
        list_pairwises_2sp = []
        with open(pairwiseFile, "r") as file:
            while (1):  # Ugly !
                name, sequence, name2, sequence2 = file.readline(), file.readline(), file.readline(), file.readline()
                if not name: break # Use assert ?
                # One locus every two lines (one pairwise couple = 4 lines) : header + sequence
                locus1 = Locus(name, sequence)
                locus2 = Locus(name2, sequence2)
                group = set([])
                group.add(locus1)
                group.add(locus2)
                list_pairwises_2sp.append(group)
        return (list_pairwises_2sp)

    # Function
    list_pairwises_allsp = []
    for file in listPairwiseFiles:
        listPairwises = getPairwiseCouple(file)
        for pairwise in listPairwises:
            list_pairwises_allsp.append(pairwise)  # all pairwises in the same 1D list
    return list_pairwises_allsp

""" Proceeds to create orthogroups by putting together pairwise couples sharing a locus.
    Iterates over the orthogroups list and tag to 'True' the pairwise couple already gathered in a group to avoid
    redondancy. Writes each orthogroup in a fasta file
    Returns an integer (a list length) """
def makeOrthogroups(list_pairwises_allsp, minspec, nb_rbh, verbose, paralogs):
    # Sub-funtions

    """ Check if a locus/group has already been treated in makeOrthogroups()
        Returns a boolean """
    def checkIfTagged(pair):
        tag = True
        for element in pair:
            if not element.getTag() and tag: # use a list comprehension maybe ?
                tag = False
        return tag

    """ True means a locus/group has already been treated in makeOrthogroups()
        A stronger code would be to implement a method inside the class Locus """
    def tagGroup(pair):
        for element in pair:
            element.tagged = True

    """ Write an orthogroup in a file """
    def writeOutputFile(orthogroup, number, naming):
        name = ""
        if naming:
            name = "orthogroup_{}_with_{}_sequences_withParalogs.fasta".format(number, len(orthogroup))
        else :
            name = "orthogroup_{}_with_{}_sequences.fasta".format(number, len(orthogroup))
        result = open(name, "w")
        with result:
            for locus in orthogroup:
                if locus.getHeader()[-1] == "\n":
                    result.write("%s" % locus.getHeader())  # write geneID
                else :
                    result.write("%s\n" % locus.Header())  # write geneID
                if locus.getSequence()[-1] == "\n":
                    result.write("%s" % locus.getSequence())  # write sequence
                else :
                    result.write("%s\n" % locus.getSequence())  # write sequence
        if naming:
            os.system("mv {} outputs_withParalogs/".format(name))
        else :
            os.system("mv {} outputs/".format(name))

    """ Builds a 2D array for a summary
        Returns a numpy 2D array """
    def countings(listOrthogroups, nb_rbh):

        """ Compute the number of species from the number of RBH comparisons
            Returns an integer """
        def compute_nbspec(nb_rbh):
    
            def sum_factorielle(x):
                n = 1
                s = 0
                while n <= x:
                    s += n
                    n += 1
                return s
            
            x = 2    
            nb_specs = 0
            while x*x - sum_factorielle(x) < nb_rbh:
                x += 1 
            return x
        #listOrthogroups.sort().reverse()
        #nblines = len(listOrthogroups[0])
        nblines = 0    
        for group in listOrthogroups:
            if len(group) > nblines:
                nblines = len(group)
        matrix = np.array([[0]*compute_nbspec(nb_rbh)]*nblines)
        # empty lines are avoided : first line of the frame is the line for minimal number of sequences in a group (>=mini)
        # for now, this feature diseappear when using numpy arrays and pandas :/

        for group in listOrthogroups:
            listSpecs = []
            for loci in group:
                if loci.getHeader()[1:3] not in listSpecs:
                    listSpecs.append(loci.getHeader()[1:3])
            matrix[len(group)-1][len(listSpecs)-1] += 1

        return matrix

    """ numpy 2D array in a nice dataframe 
        Returns a pandas 2D dataframe """
    def asFrame(matrix) :
        index = [0]*len(matrix)
        colnames = [0]*len(matrix[0])
        index = [str(i+1)+" seqs" for i in range(len(matrix))]
        colnames = [str(i+1)+" sps" for i in range(len(matrix[0]))]
        df = pd.DataFrame(matrix, index=index, columns=colnames)
        return df # Mettre une selection pour ne renvoyer que les lignes et les colonnes qui somment > 0
        #return df.loc['4 seqs':'9 seqs'].loc[:,colnames[3:]]

    # Function
    list_orthogroups = [] # That will be the orthogroup list with paralogs
    i = 1  # i is for the number of the orthogroup in its output filename
    for ortho_pair1 in list_pairwises_allsp[0:-1]:
        
        if not checkIfTagged(ortho_pair1):
            orthogroup = ortho_pair1 # the orthogroup grows as we go throught the second loop

            # check for common locus between two groups - '|' means set1.union(set2).
            for ortho_pair2 in list_pairwises_allsp[list_pairwises_allsp.index(ortho_pair1) + 1:]:
                if len(orthogroup.intersection(ortho_pair2)) != 0 and not checkIfTagged(ortho_pair2): #and not checkIfTagged(ortho_pair2) - necessary ?
                    orthogroup.update(orthogroup | ortho_pair2)
                    tagGroup(ortho_pair2)

            # Second and same loop to catch misplaced intersection
            for ortho_pair2 in list_pairwises_allsp[list_pairwises_allsp.index(ortho_pair1) + 1:]:
                if len(orthogroup.intersection(ortho_pair2)) != 0 and not checkIfTagged(ortho_pair2): #and not checkIfTagged(ortho_pair2)
                    orthogroup.update(orthogroup | ortho_pair2)
                    tagGroup(ortho_pair2)

            if len(orthogroup) > 2:
                tagGroup(ortho_pair1) # potentially useless ?                
                if len(orthogroup) >= minspec:
                    list_orthogroups.append(orthogroup)
                    if paralogs:
                        # Write orthogroup txt file with paralogs
                        writeOutputFile(orthogroup, i, True)
                    i += 1                    

    # Paralogs filtering (keeps the first encountered)
    # Quite messy ...
    # I could try to implement the paralogs filtering in the first loop to reduce computation
    list_orthogroups_format = []
    j = 1
    for nofilter_group in list_orthogroups:
        new_group = []
        species = {}
        for loci in nofilter_group:
            species[loci.getHeader()[1:3]] = False
        for loci in nofilter_group:
            if not species[loci.getHeader()[1:3]]:
                new_group.append(loci)
                species[loci.getHeader()[1:3]] = True

        if len(new_group) >= minspec: # Drop too small orthogroups            
            list_orthogroups_format.append(new_group)
            writeOutputFile(new_group, j, False)
            j += 1    
    
    # print summary table. I could also try to implement a more complex code which build and fill the frame
    # as the same time as the orthogroups are build, which would avoid to parse several times the groups list
    if verbose :
        frame = countings(list_orthogroups, nb_rbh)
        df = asFrame(frame)
        print "  Summary before paralogous filtering : \n"
        print df, "\n"

    frame2 = countings(list_orthogroups_format, nb_rbh)    
    df2 = asFrame(frame2)
    print "  Summary after paralogous filtering : \n"
    print df2

    #return len(list_orthogroups)
    return len(list_orthogroups_format)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("files", help="Input files separated by commas. Each file contains all the reciprocical best hits between a pair of species")
    parser.add_argument("minspec", help="Only keep Orthogroups with at least this number of species", type=int)
    parser.add_argument("-v", "--verbose", action="store_true", help="A supplemental summary table of orthogroups before paralogs filtering will be returned")
    parser.add_argument("-p", "--paralogs", action="store_true", help="Proceeds to write orthogroups also before paralogous filtering")
    args = parser.parse_args()

    print "*** pogs.py ***"
    print "\nBuilding of orthogroups based on pairs of genes obtained by pairwise comparisons between pairs of species."
    print "Genes are gathered in orthogroups based on the principle of transitivity between genes pairs."

    os.system("mkdir outputs")
    if args.paralogs: os.system("mkdir outputs_withParalogs")
    infiles = args.files
    listPairwiseFiles = str.split(infiles, ",")
    print "\nParsing input files ..."
    list_Locus = getListPairwiseAll(listPairwiseFiles)
    print "Creating Orthogroups ..."
    print "\n"
    nb_orthogroups = makeOrthogroups(list_Locus, args.minspec, len(listPairwiseFiles), args.verbose, args.paralogs)
    print "\n{} orthogroups have been infered from {} pairwise comparisons by RBH\n".format(nb_orthogroups, len(listPairwiseFiles))

if __name__ == "__main__":
    main()
