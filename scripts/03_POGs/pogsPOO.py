#!/usr/bin/env python
#coding: utf8
# September 2017 - Author : Victor Mataigne (Station Biologique de Roscoff - ABiMS)
# Command line : ./pogsPOO.py <list_of_input_files_separated_by_commas>

"""
Improvments to do :
    - Make the writing of orthogroups files as the same time they are made (writingOutputfiles function inside makeOrthogroups function)
        - Paralogous filtering /removal is still to do
        - Removal of sub-groups does not delete all subgroups
    - Maybe make comparaison between the final orthogroup and the current pairwise instead between the two current pairwises (that would make one loop instead of two)
"""

import sys

class Locus :
    """ Definition of a locus : header + sequence """
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence
        #self.ID = ID

    def __str__(self):
        return "{}{}".format(self.header, self.sequence)

    def __eq__(self, other):
        return self.header == other.header

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

def getPairwiseCouple(pairwiseFile):
    """ Reads an output file from the 'Pairwise' tool (AdaptSearch suite) and returns its content into a list """  
    list_pairwises_2sp = []    
    with open(pairwiseFile, "r") as file:
        while(1) : # Ugly code !
            name, sequence, name2, sequence2 = file.readline(), file.readline(), file.readline(), file.readline()
            if not name : break            
            # One locus every two lines (one pairwise couple = 4 lines) : header + sequence
            locus1 = Locus(name, sequence)   
            locus2 = Locus(name2, sequence2)        
            group = set([])
            group.add(locus1)
            group.add(locus2)            
            list_pairwises_2sp.append(group)            
    return(list_pairwises_2sp)

def getListPairwiseAll(listPairwiseFiles):
    """ Applies the getPairwiseCouple() function to a list of files and return a big list with ALL pairwises couples """
    list_pairwises_allsp = []
    for file in listPairwiseFiles:
        listPairwises = getPairwiseCouple(file)
        for pairwise in listPairwises:
            list_pairwises_allsp.append(pairwise) # all pairwises in the same list            
    return list_pairwises_allsp

def makeOrthogroups(list_pairwises_allsp):
    """ Proceeds to create orthogroups by putting together pairwise couples sharing a locus """
    list_orthogroups = []
    for ortho_pair1 in list_pairwises_allsp[0:-1]:
        orthogroup = set([])
        for ortho_pair2 in list_pairwises_allsp[list_pairwises_allsp.index(ortho_pair1)+1:]:       
            if len(ortho_pair1.intersection(ortho_pair2)) != 0 : # check for common locus between two groups. Maybe not mandatory but makes the script much faster
                orthogroup.update(ortho_pair1 | ortho_pair2) # '|' means set1.union(set2). Using set unions prevents to have duplicated loci in a group.
        if len(orthogroup) != 0 :
            list_orthogroups.append(orthogroup)
    return list_orthogroups
"""
def makeOrthogroupsRecursive(list_pairwises_allsp):
    # A recursive version to avoid removeSubGroups post-processing. Actually slower than the other one because many more iterations
    # Moreover; no stopping condition yet
    list_orthogroups = []
    len1 = len(list_pairwises_allsp)
    print len1
    for ortho_pair1 in list_pairwises_allsp[0:-1]:
        orthogroup = set([])
        for ortho_pair2 in list_pairwises_allsp[list_pairwises_allsp.index(ortho_pair1)+1:]:       
            if len(ortho_pair1.intersection(ortho_pair2)) != 0 : # check for common locus between two groups. Maybe not mandatory but makes the script much faster
                orthogroup.update(ortho_pair1 | ortho_pair2) # '|' means set1.union(set2). Using set unions prevents to have duplicated loci in a group.
        if len(orthogroup) != 0 :
            list_orthogroups.append(orthogroup)
    list_orthogroups.sort(key=len, reverse=True)
    len2 = len(list_orthogroups)
    print len2
    print "-----------"
    # Recursive
    while len2 != len1: # stop condition
        list_orthogroups = makeOrthogroupsRecursive(list_orthogroups)
    return list_orthogroups
"""

def removeSubGroups(listOfGroups):
    # Parse the list of Orthogroups to remove the subGroups (created because of an non-optimal makeOrthogroups function)
    # Not fully functional yet
    list_orthogroups = []
    listOfGroups.sort(key=len, reverse=True) # Order the list of sets by length (longer first)    
    for group1 in listOfGroups[0:-1]:
        for group2 in listOfGroups[listOfGroups.index(group1)+1:]:      
            if group1.issubset(group2) and group2 not in list_orthogroups:                
                list_orthogroups.append(group2)
            elif group2.issubset(group1) and group1 not in list_orthogroups:                
                list_orthogroups.append(group1)
                
    return list_orthogroups

def writeOutputFiles(listOfGroups):
    # Write orthogroups after the call of makeOrthogroups()
    i = 1
    for group in listOfGroups :
        name = "orthogroup_%d_with_%d_sequences.fasta" %(i, len(group))
        i += 1
        result = open(name, "w")
        with result:
            for locus in group:
                result.write("%s" %locus.getHeader()) # write geneID
                result.write("%s" %locus.getSequence()) # write sequence
    return i

def main():
    infiles = sys.argv[1]
    listPairwiseFiles = str.split(infiles, ",")
    print "Parsing input files ..."
    list_Locus = getListPairwiseAll(listPairwiseFiles)
    print "Creating Orthogroups ..."
    list_orthogroups = makeOrthogroups(list_Locus)
    #list_orthogroups = makeOrthogroupsRecursive(list_Locus)

    print "Removing unwanted sub-groups ..."
    list_orthogroups_good = removeSubGroups(list_orthogroups)
    print "%d duplicated subgroups have been removed from the orthogroups list (%d orthogroups remaining)" %(len(list_orthogroups)-len(list_orthogroups_good), len(list_orthogroups_good))
    print "writing orthogroups to fasta files ..."
    total = writeOutputFiles(list_orthogroups_good)
    print "%d orthogroups have been infered from %d pairwise comparisons" %(total-1, len(listPairwiseFiles))

if __name__ == "__main__":
    main()