#!/usr/bin/env python
#coding: utf8
# September 2017 - Author : Victor Mataigne (Station Biologique de Roscoff - ABiMS)
# Command line : ./pogsPOO.py <list_of_input_files_separated_by_commas>

"""
Improvments to do :
    - Make the writing of orthogroups files as the same time they are made (writingOutputfiles function inside makeOrthogroups function)
    - Comparisons between pairwise in makeOrthogroups function
        - maybe make comparaison between the final orthogroup and the current pairwise (that would make one loop instead of two)
        - find a way to filter duplicates at once and not after (use sets)
    - Paralogous filtering is still to do
    - Removal of sub-groups made due to the way the loop is done
    - Split classes : a class pairwise and a class orthogroup ?
"""

import sys, string

class Locus :
    """ Definition of a locus : header + sequence """
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence
        #self.ID = ID

    def __str__(self):
        return "{}{}".format(self.header, self.sequence)

    def __eq__(self, other):
        return self.header == other.header and self.sequence == other.sequence

    def __hash__(self):
        return hash((self.header, self.sequence)) #renvoie le hash du tuple

    def getHeader(self):
        return self.header

    def getSequence(self):
        return self.sequence

    """
    def __getattr__(self, header):
        return self.header

    def __getattr__(self, sequence):
        return self.sequence
    """

class Orthogroup :
    """ Gather a set of loci together. The locusList does not yet forbid items other than type 'Locus' (class Locus)"""
    def __init__(self):
        self.locusList = set()

    def __str__(self):        
        return "{}".format(self.locusList)

    def addLocus(self, locus):
        # Add a single locus to the orthogroup, even if it's already in.  
        self.locusList.add(locus)

    def addLoci(self, listLoci):
        # Add a list of loci to the orthogroup, even if it's already in.
        for loci in listLoci:
            self.locusList.add(loci)

    def removeLocus(self, locus):
        # Remove a locus from the orthogroup
        self.locusList.remove(locus)

    def getGroupLength(self):
        return len(self.locusList)

    def printList(self):
        # For each locus within the orthogroup, print its header and sequence
        for locus in self.locusList:
            print(locus)  

    def getList(self):
        return self.locusList

def getPairs(pairwise_file):
    """ Reads an output file from the 'Pairwise' tool (AdaptSearch suite) and returns its content into a list """  
    list_pairwises = []    
    with open(pairwise_file, "r") as file:
        while(1) : # Ugly code !
            name, sequence, name2, sequence2 = file.readline(), file.readline(), file.readline(), file.readline()
            if not name : break            
            # One locus every two lines (one pairwise couple = 4 lines) : header + sequence
            locus1 = Locus(name, sequence)   
            locus2 = Locus(name2, sequence2)        
            group = Orthogroup()
            group.addLocus(locus1)
            group.addLocus(locus2)            
            list_pairwises.append(group)            
    return(list_pairwises)

def getListLocus(listFiles):
    """ Applies the getPairs() function to a list of files and return a big list with ALL pairwises couples """
    listOfLocus = []
    for file in listFiles:        
        listPairwises = getPairs(file)
        for pairwise in listPairwises:
            listOfLocus.append(pairwise) # all pairwises in the same list            
    return listOfLocus

def makeOrthogroups(listOfLocus):
    """ Proceeds to create orthogroups by putting together pairwise couples sharing a locus """
    list_orthogroups = []
    
    for ortho_pair1 in listOfLocus:        
        orthogroup = Orthogroup()
        for ortho_pair2 in listOfLocus[listOfLocus.index(ortho_pair1)+1:]:
            if len(ortho_pair1.getList() & ortho_pair2.getList()) != 0 : # '&' means set1.intersection(set2)                 
                orthogroup.addLoci(ortho_pair1.getList() | ortho_pair2.getList()) # '|' means set1.union(set2)
        if orthogroup.getGroupLength() != 0 :            
            list_orthogroups.append(orthogroup)            
    return list_orthogroups

def writeOutputFiles(listOfGroups):
    # Write orthogroups after the call of makeOrthogroups() ; useless now since writing is made as the same time as the orthogroup constitution (writeAGroup function)
    i = 1
    for group in listOfGroups :        
        name = "orthogroup_%d_with_%d_sequences.fasta" %(i, group.getGroupLength())
        i += 1        
        result = open(name, "w")
        with result:
            for locus in group.getList():
                result.write("%s" %locus.getHeader()) # write geneID
                result.write("%s" %locus.getSequence()) # write sequence
    return i

def main():
    infiles = sys.argv[1]
    L1 = str.split(infiles, ",")
    print "Parsing input files ..."
    list_Locus = getListLocus(L1)    
    print "Creating Orthogroups ..."
    list_orthogroups = makeOrthogroups(list_Locus)    
    print "writing orthogroups to fasta files ..."
    total = writeOutputFiles(list_orthogroups)
    print "%d orthogroups have been infered from %d pairwise comparisons" %(total, len(L1))

if __name__ == "__main__":
    main()