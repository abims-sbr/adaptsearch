#!/usr/bin/env python
# coding: utf-8
# Author : Victor Mataigne

import string, os, sys, re, random, itertools, argparse, copy, math
import pandas as pd
import numpy as np

def buildDicts(list_codons, content, dict_genetic_code, dict_aa_classif):
    """ Build dictionaries with values to 0. These dictionaries are used as starting point for each sequence counting

    Args :
        list_codons (list of str) : all codons except codons-stop
        content (int or list) : an integer (for coutings and transitions), or an empty list (for resampling)
        dict_genetic_code (dict) : the genetic code : {'aa1': [codon1, codon2,...], ...}
        dict_aa_classif (dict) : the types of the amino-acids ( {type: [aa1, aa2, ...], ...}

    Returns :
        dico_codons, dico_aa, dico_aatypes (dicts) : keys : codons/amico-acids/amico-acids types, values : 0 or []
        dico_codons_transitions, dico_aa_transitions, dico_aatypes_transitions (dicts of dicts) :
            actually, the three first dictionaries nested as values of keys codons/amico-acids/amico-acids            
    """

    # I could have make sub-routines here.
    # the copy commands are mandatory, otherwise all dictionaries will reference the same variable(s)

    dico_codons = {}
    for codon in list_codons:
        dico_codons[codon] = copy.deepcopy(content)

    dico_aa = {}
    for aa in dict_genetic_code.keys():
        dico_aa[aa] = copy.deepcopy(content)

    dico_aatypes = {}
    for aatype in dict_aa_classif.keys():
        dico_aatypes[aatype] = copy.deepcopy(content)

    dico_codons_transitions=copy.deepcopy(dico_codons)
    for key in dico_codons_transitions.keys():
        dico_codons_transitions[key]=copy.deepcopy(dico_codons)

    dico_aa_transitions = copy.deepcopy(dico_aa)
    for key in dico_aa_transitions.keys():
        dico_aa_transitions[key]=copy.deepcopy(dico_aa)

    dico_aatypes_transitions = copy.deepcopy(dico_aatypes)
    for key in dico_aatypes_transitions.keys():
        dico_aatypes_transitions[key]=copy.deepcopy(dico_aatypes)

    return dico_codons, dico_aa, dico_aatypes, dico_codons_transitions, dico_aa_transitions, dico_aatypes_transitions

def viable(seqs, pos, method):
    """ Compute if, among a set of sequences, either at least one of the codons at the specified position has not a "-",
    or not any codon has a "-"

    Args :   
        seqs : a list (the sequences, which must all have the same size)
        pos : an integer (the positions, <= len(seqs) -3)

    Returns:
        bool
    """
    codons = [seq[pos:pos+3] for seq in seqs]
    if method == "all":
        return not all("-" in codon for codon in codons)
    elif method == "any":
        return not any("-" in codon for codon in codons)

# # # Function for codons, aa, aatypes Countings -------------------------------------------------------------------------------

def computeAllCountingsAndFreqs(seq, list_codons, init_dict_codons, init_dict_aa, init_dict_classif, dict_genetic_code, dict_aa_classif):
    """ Call all functions dedicated to the computation of occurences and frequencies of codons, amino-acids, amino-acids types

    Args : see sub-routines

    Returns : 6 dictionaries (occurences and frequencies for codons, amino-acids, amino-acids types). See sub-routines for details.
    """

    # ------ Sub-routines ------ #

    def codonsCountings(seq, list_codons, init_dict_codons):
        """ Count occurences of each codon in a sequence
        First reading frame only : input sequence is supposed to be an ORF

        Args :
            seq (str) : the sequence
            list_codons (list of str) : all codons except codons-stop
            init_dict_codons (dict) : {codon1 : 0, codon2: 0, ...}

        Return :
            codon (dict) : codons (keys) and their occurences (values) in the sequence        
        """ 

        codons = copy.deepcopy(init_dict_codons)

        l = len(seq)
        
        if l%3 == 0: max_indice = l-3
        if l%3 == 1: max_indice = l-4
        if l%3 == 2: max_indice = l-5

        for codon in range(0,max_indice+1,3):
            if "-" not in seq[codon:codon+3] :    
                codons[seq[codon:codon+3]] += 1

        return codons

    def codonsFreqs(dict_codons_counts):
        """ Computes frequencies of each codon in a sequence

        Args :        
            dict_codons (dict) : the output of codonsCounting()

        Return :
            codons_freqs (dict) : codons (keys) and their frequencies (values)
        """

        codons_freqs = {}

        for key in dict_codons_counts.keys():
            freq = float(dict_codons_counts[key])/sum(dict_codons_counts.values())
            codons_freqs[key] = freq

        return codons_freqs

    def aaCountings(dict_codons, dict_genetic_code, init_dict_aa):
        """ Count occurences of each amino-acid in a sequence, based on the countings of codons (1st ORF)

        Args :
            dict_codons (dict) : the output of codonsCounting()        
            dict_genetic_code (dict) : the genetic code : {'aa1': [codon1, codon2,...], ...}
            init_dict_aa (dict) : {aa1 : 0, aa2: 0, ...}

        Return :
            dict_aa (dict) : amino-acids (keys) and their occurences (values)
        """
        dict_aa = copy.deepcopy(init_dict_aa)        

        for key in dict_codons.keys():
            for value in dict_genetic_code.values(): 
                if key in value:            
                    dict_aa[dict_genetic_code.keys()[dict_genetic_code.values().index(value)]] += dict_codons[key]

        return dict_aa

    def aaFreqs(dict_aa_counts):
        """ Computes frequencies of each amino-acid in a sequence, based on the countings of codons (1st ORF)

        Args :    
            dict_aa_counts (dict) : the output of aaCountings()

        Return :   
            dict_aa_freqs (dict) : amino-acids (keys) and their frequencies (values)
        """

        dict_aa_freqs = {}

        for key in dict_aa_counts.keys():
            freq = float(dict_aa_counts[key])/sum(dict_aa_counts.values())
            dict_aa_freqs[key] = freq

        return dict_aa_freqs

    def aatypesCountings(dict_aa, dict_aa_classif, init_dict_classif):
        """ Computes frequencies of each amino-acid type in a sequence, based on the countings of amino-acids (1st ORF)

        Args :
        - dict_aa (dict) : the output of aaCountings() 
        - dict_aa_classif (dict) : the types of the amino-acids ( {type: [aa1, aa2, ...], ...} )
        - init_dict_classif (dict) : {'polar': 0, 'apolar': 0, ...}

        Return :
            dict_aatypes (dict) : amino-acids types (keys) and their occurences (values)
        """
        dict_aatypes = copy.deepcopy(init_dict_classif)

        for key_classif in dict_aa_classif.keys():
            for key_aa in dict_aa.keys():
                if key_aa in dict_aa_classif[key_classif]:
                    dict_aatypes[key_classif] += dict_aa[key_aa]

        return dict_aatypes

    def aatypesFreqs(dict_aatypes):
        """ Computes frequencies of each amino-acid type in a sequence, based on the countings of amino-acids (1st ORF)

        Args :
            dict_aatypes (dict) : the output of aatypesCountings()

        Return :  
            dict_aatypes_freqs : amino-acids types (keys) and their frequencies (values)
        """
        dict_aatypes_freqs = {}

        for key in dict_aatypes.keys():
            freq = float(dict_aatypes[key])/sum(dict_aatypes.values())
            dict_aatypes_freqs[key] = freq

        return dict_aatypes_freqs

    # ------ The function ------ #

    codons_c = codonsCountings(seq, list_codons, init_dict_codons)
    codons_f = codonsFreqs(codons_c)
    aa_c = aaCountings(codons_c, dict_genetic_code, init_dict_aa)
    aa_f = aaFreqs(aa_c)
    aat_c = aatypesCountings(aa_c, dict_aa_classif, init_dict_classif)
    aat_f = aatypesFreqs(aat_c)

    return codons_c, codons_f, aa_c, aa_f, aat_c, aat_f 

# # # Functions for various measures (ivywrel, ekqh...) -------------------------------------------------------------------------

def computeVarious(seq, dict_aa_counts, dict_aa_types):
    """ Call al the functions for nucleic and amino-acids sequences description

    Args : See sub-routines for details.

    Returns : 6 integer or floats. See sub-routines for details.
    """

    # ------ Sub-routines ------ #

    def nbCodons(seq):
        """ Compute the number of full codons in a sequence
        Arg : seq (str): the sequence
        Return: nb_codons (int)
        """
        l = len(seq)
        if l%3 == 0: nb_codons = l/3
        if l%3 == 1: nb_codons = (l-1)/3
        if l%3 == 2: nb_codons = (l-2)/3
        return nb_codons

    def maxIndice(seq):
        """ Compute the highest indice for parsing a sequence
        Arg : seq (str): the sequence
        Return : max_indice (int)
        """
        l = len(seq)
        if l%3 == 0: max_indice = l-3
        if l%3 == 1: max_indice = l-4
        if l%3 == 2: max_indice = l-5
        return max_indice
    
    def gc12Andgc3Count(seq, nb_codons, max_indice):
        """ Compute the frequency of gc12 in a sequence
        Arg : seq (str) : the sequence
        Return : (float)
        """

        # TO IMPROVE ? : make this computation in the codonCountigns() function to avoid parsing twice the sequence ?
        # But : involves a less readable code
       
        gc12 = 0
        gc3 = 0
        
        for i in range(0, max_indice+1,3):
            if seq[i] in ["c","g"]: gc12 += 1
            if seq[i+1] in ["c","g"]: gc12 += 1
            if seq[i+2] in ["c","g"] or seq[i+2] in ["c","g"]: gc3 += 1

        return float(gc3)/nb_codons, float(gc12)/(2*nb_codons)

    def ivywrelCount(nb_codons, dict_aa_counts):
        """ Compute the sum of occurences of amino-acids IVYWREL divided by the number of codons

        Args :
            nb_codons (int) : the number of codons in the sequence
            dict_aa_counts (dict) : the output of aaCountings()

        return : (float)
        """

        IVYWREL = 0

        for aa in ["I","V","Y","W","R","E","L"]: # Impossible to make a simple sum, in case one the aa is not in the dict keys
            if aa in dict_aa_counts.keys():
                IVYWREL += dict_aa_counts[aa]

        return float(IVYWREL)/nb_codons

    def ekqhCount(dict_aa_counts):
        """ Compute the ratio of amino-acids EK/QH
        Arg : dict_aa_counts (dict) : the output of aaCountings()
        Return : (float)
        """
        ek = 0
        qh = 0
        
        ek = dict_aa_counts["E"] + dict_aa_counts["K"]
        qh = dict_aa_counts["Q"] + dict_aa_counts["H"]
        
        if qh != 0: 
            return float(ek)/qh
        else : return ek

    def payresdgmCount(dict_aa_counts):
        """ Compute the ratio of amino-acids PASRE/SDGM
        Arg : dict_aa_counts (dict) : the output of aaCountings()
        Return : (float)
        """
        payre = 0
        sdgm = 0
        
        for aa in ["P","A","Y","R","E"]:
            payre += dict_aa_counts[aa]
        for aa in ["S","D","G","M"]:
            sdgm += dict_aa_counts[aa]

        if sdgm != 0: 
            return float(payre)/sdgm
        else : return payre

    def purineLoad(seq, nb_codons):
        """ Compute the purine load indice of a sequence
        Args : 
            seq (str) : the sequence
            nb_codons (int) : the number of codons in the sequence

        Return  (float)
        """

        # TO IMPROVE ? : make this computation in the codonCountigns() function to avoid parsing twice the sequence ?
        # But : involves a less readable code

        g12, g3, A, c12, c3, T = 0.0,0.0,seq.count("a"),0.0,0.0,seq.count("t")

        # g3 and c3 : g and c in 3rd position of a codon
        s = ""
        for i in range(2, len(seq), 3):
            s += seq[i]
        g3 = s.count("g")
        c3 = s.count("c")

        # g12 and c12 : g and c in 1st and 2d positions of a codon
        s = ""
        for i in range(0, len(seq), 3):
            s += seq[i]
        g12 = s.count("g")
        c12 = s.count("c")
        s = ""
        for i in range(1, len(seq), 3):
            s += seq[i]
        g12 += s.count("g")
        c12 += s.count("c")
        
        return float(1000*(g12+g3+A-c12-c3-T))/(3*nb_codons)

    def cvp(dict_aatypes):
        """ Compute the difference nb_charged_aamino_acids - nb_polar_amino_acids
        Return: (int)
        """
        return dict_aatypes["charged"] - dict_aatypes["polar"]

    # ------ The function ------ #
    
    nb_codons = nbCodons(seq)
    max_indice = maxIndice(seq)
    GC3, GC12 = gc12Andgc3Count(seq, nb_codons, max_indice)
    IVYWREL, EKQH, PAYRESDGM = ivywrelCount(nb_codons, dict_aa_counts), ekqhCount(dict_aa_counts), payresdgmCount(dict_aa_counts)
    purineload, CvP = purineLoad(seq, nb_codons), cvp(dict_aa_types)
    return GC3, GC12, IVYWREL, EKQH, PAYRESDGM, purineload, CvP

# # # Function for codons, aa, aatypes Transitions -----------------------------------------------------------------------------

def computeAllBiases(seq1, seq2, dico_codons_transi, dico_aa_transi, dico_aatypes_transi, reversecode, reverseclassif) :
    """ Compute all biases (transisitions codon->codon, aa->-aa, aa_type->aa_type) between two sequences    

    Args : See sub-routines for details

    Returns 3 dictionaries of dictionaries. See sub-routines for details
    """

    # ------ Sub-routines ------ #

    def codons_transitions(seq1, seq2, dico_codons_transi):
        """ Compute the number of transitions from a codon of a sequence to the codon of a second sequence

        Args :
            seq1 (str) : the first sequence
            seq2 (str) : the second sequence
            dico_codons_transi (dict of dicts) : { codon1 : {codon1: 0, codon2 : 0, ...}, ..}

        Return :
            codons_transi (dict of dicts) : the occurences of each codon to codon transition
        """

        codons_transi = copy.deepcopy(dico_codons_transi)

        for i in range(0, len(seq1), 3):
            # check if no indel and if len seq[i:i+3] is really 3 (check for the last codon)
            if viable([seq1, seq2], i, "any") and len(seq1[i:i+3]) == 3 and len(seq2[i:i+3]) == 3 :
                codons_transi[seq1[i:i+3]][seq2[i:i+3]] += 1

        return codons_transi

    def codons_transitions_freqs(codons_transitions_counts):
        """ Computes frequencies of codon transitions between two sequences

        Arg : 
            codons_transitions_counts (dict) : the output of codons_transitions()

        Return : 
            codons_transi_freqs (dict of dicts) : the frequencies of each codon to codon transition
        """

        codons_transi_freqs = {}

        for key in codons_transitions_counts.keys():
            codons_transi_freqs[key] = {}
            for key2 in codons_transitions_counts[key].keys():
                if sum(codons_transitions_counts[key].values()) != 0:
                    freq = float(codons_transitions_counts[key][key2])/sum(codons_transitions_counts[key].values())
                    codons_transi_freqs[key][key2] = freq
                else : 
                    codons_transi_freqs[key][key2] = 0 
        return codons_transi_freqs

    def aa_transitions(dico_codons_transi, dico_aa_transi, reversecode):
        """ Compute the number of transitions from an amino-acid of a sequence to the amino-acid of a second sequence

        Args :
            dico_codons_transi (dict of dicts) : the codons transitions computed by codons_transitions()
            dico_aa_transi (dict of dicts) : { aa1 : {aa1: 0, aa2 : 0, ...}, ..}
            reversecode (dict) : the reversed genetic code {aa1 : [codons],...} -> {codon1: aa1, codon2: aa2, ...}

        Return :
            aa_transi (dict of dicts) : the occurences of each aa to aa transition
        """

        aa_transi = copy.deepcopy(dico_aa_transi)
        
        for k in dico_codons_transi.keys():
            newk = reversecode[k]
            for k2 in dico_codons_transi[k].keys():
                newk2 = reversecode[k2]               
                aa_transi[newk][newk2] += dico_codons_transi[k][k2]
        
        return aa_transi

    def aa_transitions_freqs(aa_transitions_counts):
        """ Computes frequencies of amico-acids transitions between two sequences        
        Arg : aa_transitions_counts (dict of dicts): the output of aa_transitions()
        Return : aa_transi_freqs (dict of dicts) : the frequencies of each aa to aa transition
        """

        aa_transi_freqs = {}

        for key in aa_transitions_counts.keys():
            aa_transi_freqs[key] = {}
            for key2 in aa_transitions_counts[key].keys():
                if sum(aa_transitions_counts[key].values()) != 0:
                    freq = float(aa_transitions_counts[key][key2])/sum(aa_transitions_counts[key].values())
                    aa_transi_freqs[key][key2] = freq
                else :
                    aa_transi_freqs[key][key2] = 0
        return aa_transi_freqs

    def aatypes_transitions(dico_aa_transi, dico_aatypes_transi, reverseclassif):
        """ Compute the number of transitions from an amino-acid type of a sequence to the amino-acid type of a second sequence

        Args :
            dico_aa_transi (dict of dicts) : the output of aa_transitions()
            dico_aatypes_transi (dict of dicts) : { type1 : {type1: 0, type2 : 0, ...}, ..}
            reverseclassif (dict) : the reversed amino-acid clasification {aa1: type, aa2: type, ...}

        Return :
            aatypes_transi (dict of dicts) : the occurences of each aatype to aatype transition
        """

        aatypes_transi = copy.deepcopy(dico_aatypes_transi)
        for k in dico_aa_transi.keys():
            newk = reverseclassif[k]
            for k2 in dico_aa_transi[k].keys():
                newk2 = reverseclassif[k2]
                aatypes_transi[newk][newk2] += dico_aa_transi[k][k2]

        return aatypes_transi 

    def aatypes_transitions_freqs(aatypes_transitions_counts):
        """ Computes frequencies of amico-acids types transitions between two sequences
        Args : aatypes_transitions_counts (dict of dicts) : the output of aatypes_transitions()
        Return : aatypes_transi_freqs (dict of dicts) : the frequencies of each aatype to aatype transition
        """

        aatypes_transi_freqs = {}

        for key in aatypes_transitions_counts.keys():
            aatypes_transi_freqs[key] = {}
            for key2 in aatypes_transitions_counts[key].keys():
                if sum(aatypes_transitions_counts[key].values()) != 0:
                    freq = float(aatypes_transitions_counts[key][key2])/sum(aatypes_transitions_counts[key].values())
                    aatypes_transi_freqs[key][key2] = freq
                else :
                    aatypes_transi_freqs[key][key2] = 0
        return aatypes_transi_freqs

    


    # ------ The function ------ #

    codons_transitions = codons_transitions(seq1, seq2, dico_codons_transi)
    codons_transitions_freqs = codons_transitions_freqs(codons_transitions)
    aa_transitions = aa_transitions(codons_transitions, dico_aa_transi, reversecode)
    aa_transitions_freqs = aa_transitions_freqs(aa_transitions)
    aatypes_transitions = aatypes_transitions(aa_transitions, dico_aatypes_transi, reverseclassif)
    aatypes_transitions_freqs = aatypes_transitions_freqs(aatypes_transitions)
    
    return codons_transitions, codons_transitions_freqs, aa_transitions, aa_transitions_freqs, aatypes_transitions, aatypes_transitions_freqs

def all_sed(codons_c, aa_c, aat_c, codons_transitions, aa_transitions, aatypes_transitions, dico_codons_transi, dico_aa_transi, dico_aatypes_transi):

    def compute_sed(transi, counts, dico):
        """ Compute the substitution exchangeability disequilibrium (SED) from one species A to another B between codons/aa//aatypes couples

        Args:
            transi ; dict - dictionaries of all counts of transition from codon/aa/aatype X to Y from sp A to sp B
            counts : dict - dictionaries of codons/aa/aatypes counts in species A
            dico : dict - a dictionary (nested) with all values to 0

        """
        dict_sed = copy.deepcopy(dico)

        for key in transi.keys():
            for key2 in transi.keys():                
                if counts[key] != 0 and float(transi[key2][key])/counts[key] != 0.0:
                    x = (float(transi[key][key2])/counts[key]) / (float(transi[key2][key])/counts[key])
                    dict_sed[key][key2] = - pow(2,1-x)+1
                else :
                    dict_sed[key][key2] = 'NA'

        return dict_sed

    codons_sed = compute_sed(codons_transitions, codons_c, dico_codons_transi)
    aa_sed = compute_sed(aa_transitions, aa_c, dico_aa_transi)
    aatypes_sed = compute_sed(aatypes_transitions, aat_c, dico_aatypes_transi)

    return codons_sed, aa_sed, aatypes_sed

# # # Function for random resampling --------------------------------------------------------------------------------------------

def sampling (dict_seq, nb_iter, len_sample, list_codons, genetic_code, aa_classif, reversecode, reverseclassif):
    """ Resample randomly codons from sequences (sort of bootsrap)

    Args :
        dict_seq (dict) : contains the species name and sequences used ( { 'sp1': seq1, 'sp2': seq2, ...}) without '-' removal
        nb_iter (int) : number of resampling iterations (better >= 1000)
        len_sample (int) : length (in codons) of a resampled sequence (better >= 1000) 
        list_codons (list of str): all codons except codons-stop
        genetic_code (dict) : the genetic code : {'aa1': [codon1, codon2,...], ...}
        aa_classif (dict) : the types of the amino-acids : {type: [aa1, aa2, ...], ...}
        reversecode (dict) : the reversed genetic code : {codon1: aa1, codon2: aa2, ...}
        reverseclassif (dict) : the reversed amino-acid : clasification {aa1: type, aa2: type, ...}

    Returns :
        codons_lst, aa_lst, classif_lst (dicts) : keys : codons/aa/aatypes, values : []        
        codons_transitions_lst, aa_transitions_lst, classif_transitions_lst (dict of dicts) : keys : codons/aa/aatypes, values : the 3 previous dicts
    """

    # Initialize empty dictionaries for countings and transitions. It's also possible to isntanciate these ones in the main() but it would make a function with ~14 parameters
    codons_0, aa_0, classif_0, codons_transitions_0, aa_transitions_0, classif_transitions_0 = buildDicts(list_codons, 0, genetic_code, aa_classif)
    codons_lst, aa_lst, classif_lst, codons_transitions_lst, aa_transitions_lst, classif_transitions_lst = buildDicts(list_codons, [], genetic_code, aa_classif)
    
    # determine the max position where sampling is possible
    l = len(dict_seq.values()[1])
    if l%3 == 0: max_indice = l-3
    if l%3 == 1: max_indice = l-4
    if l%3 == 2: max_indice = l-5

    # List of positions to resample
    viable_positions = [pos for pos in range(0,max_indice,3) if viable(dict_seq.values(), pos, "all")]
    sample_positions = np.random.choice(viable_positions, len_sample)

    # nb_iter resampled sequences
    for i in range(nb_iter):
        if (i+1)%(nb_iter/10) == 0:
            print "    "+str( (i+1)*100/nb_iter)+"%"

        seqa, seqb = "", ""
        for pos in sample_positions:
            codona, codonb = "---", "---"
            # The sequence to be resampled in this position is randomly chosen  ; no "-" resampled
            while "-" in codona :
                codona = dict_seq.values()[random.randrange(0, len(dict_seq.keys())-1)][pos:pos+3]
            while "-" in codonb :
                codonb = dict_seq.values()[random.randrange(0, len(dict_seq.keys())-1)][pos:pos+3]
            seqa += codona
            seqb += codonb
        
        # dictionaries : frequences of codons, aa, aatypes (seq1)
        codons_occ_tmp, codons_freq_tmp, aa_occ_tmp, aa_freq_tmp, aatypes_occ_tmp, aatypes_freq_tmp = computeAllCountingsAndFreqs(seqa, list_codons, codons_0, aa_0, classif_0, genetic_code, aa_classif)
        # dictionaries frequences of transitions (seqa->seqb)
        codons_transitions_tmp, codons_transitions_freq_tmp, aa_transition_tmp, aa_transitions_freq_tmp, aatypes_transitions_tmp, aatypes_transitions_freq_tmp = computeAllBiases(seqa, seqb, codons_transitions_0, aa_transitions_0, classif_transitions_0, reversecode, reverseclassif)
        
        # Adding occurrences in final dicts
        for key in codons_freq_tmp.keys():
            codons_lst[key].append(codons_freq_tmp[key])
        for key in aa_freq_tmp.keys():
            aa_lst[key].append(aa_freq_tmp[key])
        for key in aatypes_freq_tmp.keys():
            classif_lst[key].append(aatypes_freq_tmp[key])
        
        # Adding occurrences in final dicts (transitions)
        for key in codons_transitions_freq_tmp.keys():
            for key2 in codons_transitions_freq_tmp[key].keys():
                codons_transitions_lst[key][key2].append(codons_transitions_freq_tmp[key][key2])
        for key in aa_transitions_freq_tmp.keys():
            for key2 in aa_transitions_freq_tmp[key].keys():
                aa_transitions_lst[key][key2].append(aa_transitions_freq_tmp[key][key2])
        for key in aatypes_transitions_freq_tmp.keys():
            for key2 in aatypes_transitions_freq_tmp[key].keys():
                classif_transitions_lst[key][key2].append(aatypes_transitions_freq_tmp[key][key2])
    
    return codons_lst, aa_lst, classif_lst, codons_transitions_lst, aa_transitions_lst, classif_transitions_lst

def testPvalues(dict_counts, dict_resampling, nb_iter, method):
    """ Computes where the observed value is located in the expected counting distribution

    Args :
        dict_counts (dict) : observed frequencies obtained from the functions computeAllCountingsAndFreqs() or computeAllBiases()
        dict_resampling (dict) : expected frequencies obtained from the functions computeAllCountingsAndFreqs() or computeAllBiases() within the sampling() function
    
    Return :
        pvalue (dict, dict of dicts) : the pvalues of all observed countings (dict) and transitions (dict of dicts)
    

    pnorm computes the pvalue to have a value inferior to the observed value under a normal distribution
    One sided to left tail : 
        p < 0.05 indicates significantly lower counts
        p > 0.95 indicates significantly higher counts
    """


    def p_resampling(obs, values, nb_iter):
        """ The pvalue is the proportion of bootsrapped values smaller than the observed value
        If p = 0.025 : 2.5% of the bootstrapped values are smaller than the observed value
           p < 0.025 : the obs value is most likely significantly lower.        
        If p = 0.975 : 97.5% of the bootstrapped values are smaller than the observed value
           p > 0.975   the obs value is most likely significantly higher.        

        Args :
            obs : int or float - the observed value
            values : list - values of resampling (int or floats)
            nb_iter : int - the number of resampled values (=len(values))

        Return :
            pvalue (float)
        """

        num = len([x for x in values if x < obs])
        return float(num + 1) / (nb_iter+1)

    def testPvalue(obs, exp, nb_iter):
        """ Compute a pvalue

        Args :
            exp (list of floatsà : a list of length nb_iter, containing expected frequencies of a codon/aa/aatype at each iteration
            obs (float) : the observed value
            nb_iter (int) : the number of iterations for resampling

        Returns :
            pvalue (float)
        """

        max_val = nb_iter-1
        min_val = 0
        test_val = (max_val+min_val)/2
          
        while max_val-min_val > 1:
            if obs > exp[test_val]:
                min_val = test_val
                test_val = (max_val+min_val)/2
            elif obs < exp[test_val]:
                max_val = test_val
                test_val = (max_val+min_val)/2
            else:
                break

        pvalue = float(test_val+1)/(nb_iter+1)

        return pvalue

    # ------ The function ------ #

    pvalues = {}

    for key in dict_resampling.keys():
        if type(dict_resampling.values()[1]) is not dict :
            if method == 'origin':
                pvalues[key] = testPvalue(dict_counts[key], dict_resampling[key], nb_iter)
            elif method == 'pnorm':
                pvalues[key] = scipy.stats.norm.cdf(dict_counts[key], np.mean(dict_resampling[key]), np.std(dict_resampling[key]))
            elif method == 'p_resampling':
                pvalues[key] = p_resampling(dict_counts[key], dict_resampling[key], nb_iter)
        else :
            pvalues[key] = {}
            for key2 in dict_resampling[key].keys():
                if method == 'origin':
                    pvalues[key][key2] = testPvalue(dict_counts[key][key2], dict_resampling[key][key2], nb_iter)
                elif method == 'pnorm':
                    pvalues[key][key2] = scipy.stats.norm.cdf(dict_counts[key][key2], np.mean(dict_resampling[key][key2]), np.std(dict_resampling[key][key2]))
                elif method == 'p_resampling':
                    pvalues[key][key2] = p_resampling(dict_counts[key][key2], dict_resampling[key][key2], nb_iter)

    return pvalues

def main():
    
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("sequences_file", help="File containing sequences (the output of the tool 'ConcatPhyl'")
    parser.add_argument("considered_species", help="The species name, separated by commas (must be the same than in the sequences_file). It is possible to consider only a subset of species.")
    parser.add_argument("species_for_bootstrap", help="The species which will be used for bootstrapping, separated by commas. It is possible to consider only a subset of species.")
    parser.add_argument("iteration", help="The number of iterations for bootstrapping (better if => 1000)", type=int)
    parser.add_argument("sample_length", help="The lenght of a bootstrapped sequences (better if >= 1000", type=int)
    args = parser.parse_args()

    print "\n ------ Occurences and frequencies of codons, amino-acids, amino-acids types -------\n"

    print "The script counts the number of codons, amino acids, and types of amino acids in sequences,"
    print "as well as the mutation bias from one item to another between 2 sequences.\n"

    print "Counting are then compared to empirical p-values, obtained from bootstrapped sequences obtained from a subset of sequences."
    print "In the output files, the pvalues indicate the position of the observed data in a distribution of empirical countings obtained from"
    print "a resample of the data. Values above 0.95 indicate a significantly higher counting, values under 0.05 a significantly lower counting."

    print "    Sequences file : {}".format(args.sequences_file)
    print "    Species retained for countings : {}\n".format(args.considered_species)

    print "Processing : reading input file, opening output files, building dictionaries."
    
    # make pairs
    list_species = str.split(args.considered_species, ",")
    list_species_boot = str.split(args.species_for_bootstrap, ",")
    pairs_list=list(itertools.combinations(list_species,2))    

    # read sequences
    sequences_for_counts = {}
    sequences_for_resampling = {}
    with open(args.sequences_file, "r") as file:
        for line1,line2 in itertools.izip_longest(*[file]*2):
            species = line1.strip(">\r\n")
            sequence = line2.strip("\r\n")
            if species in list_species:                
                sequences_for_counts[species] = sequence
            if species in list_species_boot:
                sequences_for_resampling[species] = sequence

    print "    Warning : countings might be biased and show high differences between species because of high variations of the indels proportions among sequences."
    print "    Frequences are more representative."

    print "\n    Indels percent :"

    for k,v in sequences_for_counts.items():        
        print "    {} : {} %".format(k, float(v.count("-"))/len(v)*100)

    # useful dictionaries
    dict_genetic_code={"F":["ttt","ttc"],
                       "L":["tta","ttg","ctt","ctc","cta","ctg"],
                       "I":["att","atc","ata"],
                       "M":["atg"],
                       "V":["gtt","gtc","gta","gtg"],
                       "S":["tct","tcc","tca","tcg","agt","agc"],
                       "P":["cct","cca","ccg","ccc"],
                       "T":["act","acc","aca","acg"],
                       "A":["gct","gcc","gca","gcg"],
                       "Y":["tat","tac"],
                       "H":["cat","cac"],
                       "Q":["caa","cag"],
                       "N":["aat","aac"],
                       "K":["aaa","aag"],
                       "D":["gat","gac"],
                       "E":["gaa","gag"],
                       "C":["tgt","tgc"],
                       "W":["tgg"],
                       "R":["cgt","cgc","cga","cgg","aga","agg"],
                       "G":["ggt","ggc","gga","ggg"]}

    dict_aa_classif={"unpolar":["G","A","V","L","M","I"],
                     "polar":["S","T","C","P","N","Q"],
                     "charged":["K","R","H","D","E"],
                     "aromatics":["F","Y","W"]}

    reversecode={v:k for k in dict_genetic_code for v in dict_genetic_code[k]}
    reverseclassif={v:k for k in dict_aa_classif for v in dict_aa_classif[k]}

    # codons list (without stop codons)
    nucleotides = ['a', 'c', 'g', 't']
    list_codons = [''.join(comb) for comb in itertools.product(nucleotides, repeat=3)]
    list_codons.remove("taa")
    list_codons.remove("tag")
    list_codons.remove("tga")
    
    # Store already computed species + row.names in output files
    index = []
    index_transi = []

    # Final dictionaries writed to csv files
    all_codons = {}
    all_aa = {}
    all_aatypes = {}
    all_various = {}
    all_codons_transitions = {} # Not used because too much : 61*61 columns
    all_aa_transitions = {}
    all_aatypes_transitions = {}
    
    # RUN

    print "\nProcessing : resampling ..."
    print "    Parameters : {niter} iterations, {lensample} codon per resampled sequence, species used : {species}\n".format(niter=args.iteration, lensample=args.sample_length, species=args.species_for_bootstrap)

    codons_boot, aa_boot, aatypes_boot, codons_transi_boot, aa_transi_boot, aatypes_transi_boot = sampling(sequences_for_resampling, args.iteration, args.sample_length, list_codons, dict_genetic_code, dict_aa_classif, reversecode, reverseclassif)
    print "    Done.\n"

    print "Processing : countings....\n"

    # Initialize empty dictionaries for countings and transitions
    init_dict_codons, init_dict_aa, init_dict_classif, dico_codons_transitions, dico_aa_transitions, dico_aatypes_transitions = buildDicts(list_codons, 0, dict_genetic_code, dict_aa_classif)

    for pair in pairs_list:
        p1, p2 = pair[0], pair[1]
        if p1 not in index:
            print "Countings on {}".format(p1)

            p1_codons_counts, p1_codons_freqs, p1_aa_counts, p1_aa_freqs, p1_aatypes_counts, p1_aatypes_freqs = computeAllCountingsAndFreqs(sequences_for_counts[p1], list_codons, init_dict_codons, init_dict_aa, init_dict_classif, dict_genetic_code, dict_aa_classif)
            p1_GC3, p1_GC12, p1_IVYWREL, p1_EKQH, p1_PAYRESDGM, p1_purineload, p1_CvP = computeVarious(sequences_for_counts[p1], p1_aa_counts, p1_aatypes_freqs)

            
            p1_codons_pvalues = testPvalues(p1_codons_freqs, codons_boot, args.iteration, 'p_resampling')
            p1_aa_pvalues = testPvalues(p1_aa_freqs, aa_boot, args.iteration, 'p_resampling')
            p1_aatypes_pvalues = testPvalues(p1_aatypes_freqs, aatypes_boot, args.iteration, 'p_resampling')

            all_codons[p1+"_obs_counts"] = p1_codons_counts
            all_codons[p1+"_obs_freqs"] = p1_codons_freqs
            all_codons[p1+"_pvalues"] = p1_codons_pvalues
            all_aa[p1+"_obs_counts"] = p1_aa_counts
            all_aa[p1+"_obs_freqs"] = p1_aa_freqs
            all_aa[p1+"_pvalues"] = p1_aa_pvalues
            all_aatypes[p1+"_obs_counts"] = p1_aatypes_counts
            all_aatypes[p1+"_obs_freqs"] = p1_aatypes_freqs
            all_aatypes[p1+"_pvalues"] = p1_aatypes_pvalues
            all_various[p1] = [p1_GC3, p1_GC12, p1_IVYWREL, p1_EKQH, p1_PAYRESDGM, p1_purineload, p1_CvP]

            index.append(p1)

        if p2 not in index:
            print "Countings on {}".format(p2)
            
            p2_codons_counts, p2_codons_freqs, p2_aa_counts, p2_aa_freqs, p2_aatypes_counts, p2_aatypes_freqs = computeAllCountingsAndFreqs(sequences_for_counts[p2], list_codons, init_dict_codons, init_dict_aa, init_dict_classif, dict_genetic_code, dict_aa_classif)
            p2_GC3, p2_GC12, p2_IVYWREL, p2_EKQH, p2_PAYRESDGM, p2_purineload, p2_CvP = computeVarious(sequences_for_counts[p2], p2_aa_counts, p2_aatypes_freqs)
            
            p2_codons_pvalues = testPvalues(p2_codons_freqs, codons_boot, args.iteration, 'p_resampling')
            p2_aa_pvalues = testPvalues(p2_aa_freqs, aa_boot, args.iteration, 'p_resampling')
            p2_aatypes_pvalues = testPvalues(p2_aatypes_freqs, aatypes_boot, args.iteration, 'p_resampling')

            all_codons[p2+"_obs_counts"] = p2_codons_counts
            all_codons[p2+"_obs_freqs"] = p2_codons_freqs
            all_codons[p2+"_pvalues"] = p2_codons_pvalues
            all_aa[p2+"_obs_counts"] = p2_aa_counts
            all_aa[p2+"_obs_freqs"] = p2_aa_freqs
            all_aa[p2+"_pvalues"] = p2_aa_pvalues
            all_aatypes[p2+"_obs_counts"] = p2_aatypes_counts
            all_aatypes[p2+"_obs_freqs"] = p2_aatypes_freqs
            all_aatypes[p2+"_pvalues"] = p2_aatypes_pvalues
            all_various[p2] = p2_GC3, p2_GC12, p2_IVYWREL, p2_EKQH, p2_PAYRESDGM, p2_purineload, p2_CvP

            index.append(p2)

        if (p1, p2) not in index_transi and p1 in sequences_for_counts and p2 in sequences_for_counts:
            print "Countings transitions between {} and {}".format(p1, p2)
            codons_transitions, codons_transitions_freqs, aa_transitions, aa_transitions_freqs, aatypes_transitions, aatypes_transitions_freqs = computeAllBiases(sequences_for_counts[p1], sequences_for_counts[p2], dico_codons_transitions, dico_aa_transitions, dico_aatypes_transitions, reversecode, reverseclassif)
            
            # Ajout
            codons_sed, aa_sed, aatypes_sed = all_sed(p1_codons_counts, p1_aa_counts, p1_aatypes_counts, codons_transitions, aa_transitions, aatypes_transitions, dico_codons_transitions, dico_aa_transitions, dico_aatypes_transitions)

            index_transi.append((p1,p2))

            p1p2_codons_pvalues = testPvalues(codons_transitions_freqs, codons_transi_boot, args.iteration, 'p_resampling')
            p1p2_aa_pvalues = testPvalues(aa_transitions_freqs, aa_transi_boot, args.iteration, 'p_resampling')
            p1p2_aatypes_pvalues = testPvalues(aatypes_transitions_freqs, aatypes_transi_boot, args.iteration, 'p_resampling')

            all_codons_transitions[p1+">"+p2+"_obs_counts"] = codons_transitions
            all_codons_transitions[p1+">"+p2+"_obs_freqs"] = codons_transitions_freqs
            all_codons_transitions[p1+">"+p2+"_pvalues"] = p1p2_codons_pvalues
            all_aa_transitions[p1+">"+p2+"_obs_counts"] = aa_transitions
            all_aa_transitions[p1+">"+p2+"_obs_freqs"] = aa_transitions_freqs
            all_aa_transitions[p1+">"+p2+"_pvalues"] = p1p2_aa_pvalues
            all_aatypes_transitions[p1+">"+p2+"_obs_counts"] = aatypes_transitions
            all_aatypes_transitions[p1+">"+p2+"_obs_freqs"] = aatypes_transitions_freqs            
            all_aatypes_transitions[p1+">"+p2+"_pvalues"] = p1p2_aatypes_pvalues

            all_codons_transitions[p1+">"+p2+"_sed"] = codons_sed
            all_aa_transitions[p1+">"+p2+"_sed"] = aa_sed
            all_aatypes_transitions[p1+">"+p2+"_sed"] = aatypes_sed

            index_transi.append((p1, p2))

    print "\n    Done.\n"

    print "Processing : creating dataframes ..."    

    frame_codons = pd.DataFrame(all_codons).T.astype('object')
    frame_aa = pd.DataFrame(all_aa).T.astype('object')    
    frame_aatypes = pd.DataFrame(all_aatypes).T.astype('object')

    frame_codons_transitions = pd.concat({k: pd.DataFrame(v) for k, v in all_codons_transitions.items()}).unstack()
    frame_codons_transitions.columns = frame_codons_transitions.columns.map('>'.join)

    frame_aa_transitions = pd.concat({k: pd.DataFrame(v) for k, v in all_aa_transitions.items()}).unstack()
    frame_aa_transitions.columns = frame_aa_transitions.columns.map('>'.join)

    frame_aatypes_transitions = pd.concat({k: pd.DataFrame(v) for k, v in all_aatypes_transitions.items()}).unstack()
    frame_aatypes_transitions.columns = frame_aatypes_transitions.columns.map('>'.join)

    frame_various = pd.DataFrame(all_various).T
    frame_various.columns = ["GC3","GC12","IVYWREL","EKQH","PAYRESDGM","purineload", "CvP"]

    frame_codons.index.name, frame_aa.index.name, frame_aatypes.index.name = "Species", "Species","Species"
    frame_aa_transitions.index.name, frame_aatypes_transitions.index.name, frame_various.index.name = "Species","Species","Species"

    print "Writing dataframes to output files ...\n"

    frame_codons.to_csv("codons_freqs.csv", sep=",", encoding="utf-8")
    frame_aa.to_csv("aa_freqs.csv", sep=",", encoding="utf-8")
    frame_aatypes.astype('object').to_csv("aatypes_freqs.csv", sep=",", encoding="utf-8")
    frame_codons_transitions.to_csv("codons_transitions_freqs.csv", sep=",", encoding="utf-8")
    frame_aa_transitions.to_csv("aa_transitions_freqs.csv", sep=",", encoding="utf-8")
    frame_aatypes_transitions.to_csv("aatypes_transitions_freqs.csv", sep=",", encoding="utf-8")
    frame_various.to_csv("gc_and_others_freqs.csv", sep=",", encoding="utf-8")

    print "Done."

if __name__ == "__main__":
    main()