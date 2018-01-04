#!/usr/bin/env python
# coding: utf8
# September 2017 - Author : Victor Mataigne (Station Biologique de Roscoff - ABiMS)
# Command line : ./pogsPOO.py <list_of_input_files_separated_by_commas> <minimal number of species per group>

"""
What it does:
    - pogs.py parses output files from the "pairwise" tool of the AdaptSearch suite and proceeds to gather genes in orthogroups (using transitivity).
    - A minimal number of species per group can be set.

Improvments to do :
    - Paralogous filtering /removal is still to do

Errors corrected (need to test a bit more to be sure)
    - makeOrthogroups is still incorrect
        - ERROR 1
            - Pour peu qu'une paire soit associee au groupe grace a une paire qui sera parcourue après elle, cette premiere paire n'est pas ajoutee
            - PAS SUR QUE CA MARCHE dans tous les cas : pb regle en faisant la boucle deux fois de suite        
        - ERROR 2         
            - Comportement bugue de intersection : 
                La longueur de l'intersect est parfois nulle alors qu'on a bien un loci1.getHeader() == loci2.getHeader() à True (j'ai teste chaque paire)
                On dirait qu'il faut créer une classe Orthogroup qui hérite de set et redéfinir la méthode intersection ...
"""

import os, argparse

class Locus:
    """ Definition of a locus : header + sequence """

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

def getListPairwiseAll(listPairwiseFiles):
    """ Applies the getPairwiseCouple() function to a list of files and return a big list with ALL pairwises couples """
    
    # Sub-Function
    def getPairwiseCouple(pairwiseFile):
        """ Reads an output file from the 'Pairwise' tool (AdaptSearch suite) and returns its content into a list """
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
            list_pairwises_allsp.append(pairwise)  # all pairwises in the same list
    return list_pairwises_allsp

def makeOrthogroups(list_pairwises_allsp, minspec, args):
    """ Proceeds to create orthogroups by putting together pairwise couples sharing a locus.
     Iterates over the orthogroups list and tag to 'True' the pairwise couple already gathered in a group to avoid
     redondancy. Writes each orthogroup in a fasta file"""

    # Sub-funtions
    def checkIfTagged(pair):
        """ Check if a locus/group has already been treated in makeOrthogroups()"""
        tag = True
        for element in pair:        
            if not element.getTag() and tag: # use a list comprehension maybe ?
                tag = False
        return tag

    def tagGroup(pair):
        """ True means a locus/group has already been treated in makeOrthogroups()
         A stronger code would be to implement a method inside the class Locus"""
        for element in pair:
            element.tagged = True

    def writeOutputFile(orthogroup, number):
        """ Write an orthogroup in a file """
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
        os.system("mv {} outputs/".format(name))

    # Function
    list_orthogroups = []
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
                    writeOutputFile(orthogroup, i)
                    i += 1
    
    # print groups in the log file
    if args.verbose > 0:
        for group in list_orthogroups:
            print "*** Orthogroup with {} sequences ***".format(len(group))
            for locus in group:
                if args.verbose == 1 :
                    locus.prettyPrint()
                elif args.verbose == 2 :
                    locus.prettyPrint2()
            print ""
    
    return len(list_orthogroups)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("files", help="Input files separated by commas. Each file contains all the reciprocical best hits between a pair of species")
    parser.add_argument("minspec", help="Only keep Orthogroups with at least this number of species", type=int)
    parser.add_argument("-v", "--verbose", type=int, choices=[0,1,2], help="Content of orthogroups will be displayed on the screen")
    args = parser.parse_args()

    print "*** pogs.py ***"
    print "\nBuilding of orthogroups based on pairs of genes obtained by pairwise comparisons between pairs of species."
    print "Genes are gathered in orthogroups based on the principle of transitivity between genes pairs."

    os.system("mkdir outputs")
    infiles = args.files
    listPairwiseFiles = str.split(infiles, ",")
    print "\nParsing input files ..."
    list_Locus = getListPairwiseAll(listPairwiseFiles)
    print "Creating Orthogroups ..."
    print "\n"
    nb_orthogroups = makeOrthogroups(list_Locus, args.minspec, args) # args must be passed as a parameter ; maybe not the best way to do it ?
    print "{} orthogroups have been infered from {} pairwise comparisons by RBH".format(nb_orthogroups, len(listPairwiseFiles))

if __name__ == "__main__":
    main()
