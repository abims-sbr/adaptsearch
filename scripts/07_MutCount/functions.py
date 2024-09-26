#!/usr/bin/env python
#coding: utf-8

import itertools, os

def dico(fasta_file, path_in):
    """
    Stores a fasta file in a dictionary : key/value -> header/sequence

    Args:
        - fasta_file (String) : the name of fasta file
        - path_in (String) : path to the fasta file

    Return:
        - bash1 (dict) : the dictionary header/sequence        
    """
    bash1 = {}    

    with open(path_in+'/'+fasta_file, 'r') as F1:
        for h,s in itertools.zip_longest(*[F1]*2):            
            fasta_name = h[1:3]
            sequence = s[:-1]
            if fasta_name not in list(bash1.keys()):
                bash1[fasta_name] = sequence
            else:
                print(fasta_name)
   
    return bash1 # same length for all (alignment)

def write_output(names, sps_list, out_dir, results_dict):
    """ Write results in csv files. There is one file per counted element (one file per amino-acid, one file per indice ...)

    Args:
        - names (list) : list with the names of elems
        - sps_list (list) : species names, sorted alphabetically
        - out_dir (String) : output directory
        - results_dict (dict) : vcounts values of each element for each input file (keys names : elems from 'names argument')

    """
    for name in names:
        out = open(name+".csv", 'w')
        out.write('Group,' + sps_list[0:-1]+'\n')
        for group in list(results_dict.keys()):
            count_of_elems = ''
            for specs in sorted(results_dict[group].keys()):
                count_of_elems += str(results_dict[group][specs][name]) + ','
            out.write(group + ',' + count_of_elems[0:-1] + '\n')
        out.close()
        os.system('mv %s.csv %s/' %(name, out_dir))

def fill_with_NaN(what):
    """ Used to create a dict only with NaN values ; used when a species is not present in an orthogroup

    Args:
        - what (list of Strings) : the names of the elements studied (nucleotide, amino-acids, indices of thermostability ...)

    Return:
        - NaN_values (dict) : dictionary with keys=elems of what, values=NaN
    """

    NaN_values = {}
    for elem in what:
        NaN_values[elem] = 'NaN'

    return NaN_values
