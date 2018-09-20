#!/usr/bin/env python
#coding: utf-8
#Author : Eric Fontanillas (2010) - Victor Mataigne (2018)

import itertools, argparse, os
import pandas as pd

## Generates bash, with key = fasta name; value = sequence (WITH GAP, IF ANY, REMOVED IN THIS FUNCTION)
def dico(fasta_file, path_in):
    '''
    Stores a fasta file in a dictionary : key/value -> header/sequence
    '''
    bash1 = {}    

    with open(path_in+'/'+fasta_file, 'r') as F1:
        for h,s in itertools.izip_longest(*[F1]*2):            
            fasta_name = h[1:3]
            sequence = s[:-1]
            if fasta_name not in bash1.keys():
                bash1[fasta_name] = sequence
            else:
                print fasta_name
   
    return(bash1, len(bash1[bash1.keys()[0]])) # same length for all (alignment)

def one_aa_counts(seq, aa):
    count = seq.count(aa)
    prop = count/len(seq)
    return count, prop

def all_aa_counts(seq):

    aa_counts = {}
    LAA =['K','R','A','F','I','L','M','V','W','N','Q','S','T','H','Y','C','D','E','P','G']

    for aa in LAA:
        aa_counts[aa] = seq.count(aa)

    return aa_counts

def all_aa_props(seq_counts):

    aa_props = {}
    for key in seq_counts.keys():
        aa_props[key] = float(seq_counts[key]) / sum(seq_counts.values())
    return aa_props

def aa_variables_counts_and_props(aa_counts):

    aa_variables_counts = {}
    aa_variables_props = {}    
    len_seq = sum(aa_counts.values()) # length of the sequence
    
    # counts of variables
    aa_variables_counts['AC'] = aa_counts['A'] + aa_counts['C']
    aa_variables_counts['APGC'] = aa_counts['A'] + aa_counts['P'] + aa_counts['G'] + aa_counts['C']
    aa_variables_counts['AVLIM'] = aa_counts['A'] + aa_counts['V'] + aa_counts['L'] + aa_counts['I'] + aa_counts['M']
    aa_variables_counts['AVLIMFYW'] = aa_variables_counts['AVLIM'] + aa_counts['F'] + aa_counts['Y'] + aa_counts['W']
    aa_variables_counts['DE'] = aa_counts['D'] + aa_counts['E']
    aa_variables_counts['DNQTSHA'] = aa_counts['D'] + aa_counts['N'] + aa_counts['Q'] + aa_counts['T'] + aa_counts['S'] + aa_counts['H'] + aa_counts['A']
    aa_variables_counts['EK'] = aa_counts['E'] + aa_counts['K']
    aa_variables_counts['ERK'] = aa_counts['E'] + aa_counts['K'] + aa_counts['K']
    aa_variables_counts['FYMINK'] = aa_counts['F'] + aa_counts['Y'] + aa_counts['M'] + aa_counts['I'] + aa_counts['N'] + aa_counts['K']
    aa_variables_counts['FYW'] = aa_counts['F'] + aa_counts['Y'] + aa_counts['W']
    aa_variables_counts['GARP'] = aa_counts['G'] + aa_counts['A'] + aa_counts['R'] + aa_counts['P']
    aa_variables_counts['IVYWREL'] = aa_counts['I'] + aa_counts['V'] + aa_counts['Y'] + aa_counts['W'] + aa_counts['R'] + aa_counts['E'] + aa_counts['L']
    aa_variables_counts['QH'] = aa_counts['Q'] + aa_counts['H']
    aa_variables_counts['RHK'] = aa_counts['R'] + aa_counts['H'] + aa_counts['K']
    aa_variables_counts['RHKDE'] = aa_counts['R'] + aa_counts['H'] + aa_counts['K'] + aa_counts['D'] + aa_counts['E']
    aa_variables_counts['STNQ'] = aa_counts['S'] + aa_counts['T'] + aa_counts['N'] + aa_counts['Q']
    aa_variables_counts['VLIM'] = aa_counts['V'] + aa_counts['L'] + aa_counts['I'] + aa_counts['M']

    # compute proportions
    for key in aa_variables_counts.keys():
        aa_variables_props[key] = float(aa_variables_counts[key]) / float(len_seq)

    if aa_variables_counts['DNQTSHA'] != 0:
        ratio_ERK_DNQTSHA = float(aa_variables_counts['ERK'])/float(aa_variables_counts['DNQTSHA'])
    else :
        ratio_ERK_DNQTSHA = -1

    if aa_variables_counts['QH'] != 0:
        ratio_EK_QH = float(aa_variables_counts['EK'])/float(aa_variables_counts['QH'])
    else :
        ratio_EK_QH = -1
    
    if aa_variables_counts['FYMINK'] != 0:
        ratio_GARP_FYMINK = float(aa_variables_counts['EK'])/float(aa_variables_counts['FYMINK'])
    else :
        ratio_GARP_FYMINK = -1
    
    if aa_variables_counts['VLIM'] != 0:
        ratio_AC_VLIM = float(aa_variables_counts['AC'])/float(aa_variables_counts['VLIM'])
        ratio_APGC_VLIM =  float(aa_variables_counts['APGC'])/float(aa_variables_counts['VLIM'])
    else :
        ratio_AC_VLIM = -1
        ratio_APGC_VLIM = -1

    aa_variables_counts['ratio_ERK_DNQTSHA'] = ratio_ERK_DNQTSHA
    aa_variables_counts['ratio_EK_QH'] = ratio_EK_QH
    aa_variables_counts['ratio_GARP_FYMINK'] = ratio_GARP_FYMINK
    aa_variables_counts['ratio_AC_VLIM'] = ratio_AC_VLIM
    aa_variables_counts['ratio_APGC_VLIM'] = ratio_APGC_VLIM

    return aa_variables_counts, aa_variables_props

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("species_list", help="List of species sperated by commas")    
    args = parser.parse_args()

    LAA =['K','R','A','F','I','L','M','V','W','N','Q','S','T','H','Y','C','D','E','P','G']
    LV = ['IVYWREL','EK','ERK','DNQTSHA','QH','ratio_ERK_DNQTSHA','ratio_EK_QH','FYMINK','GARP',
          'ratio_GARP_FYMINK','AVLIM','FYW','AVLIMFYW','STNQ','RHK','DE','RHKDE','APGC','AC',
          'VLIM','ratio_AC_VLIM','ratio_APGC_VLIM']

    path_inputs = '01_input_files'
    list_inputs = os.listdir(path_inputs)

    lsp = args.species_list.split(',')
    lsp = sorted(lsp)
    flsp = ''
    for el in lsp:
        flsp += el+','    

    path_outputs_1 = '02_tables_per_aa'
    path_outputs_2 = '02_tables_per_variable'
    os.mkdir(path_outputs_1)
    os.mkdir(path_outputs_2)

    dict_for_files_aa_counts = {}
    dict_for_files_aa_props = {}
    dict_for_files_variables_counts = {}
    dict_for_files_variables_props = {}

    # All counts and props
    for file in list_inputs:
        content, l = dico(file, path_inputs)
        aa_counts_per_seq = {}
        aa_props_per_seq = {}
        aa_variables_counts_per_seq = {}
        aa_variables_props_per_seq = {}

        for key in content.keys():
            aa_counts_per_seq[key] = all_aa_counts(content[key])
            aa_props_per_seq[key] = all_aa_props(aa_counts_per_seq[key])
            aa_variables_counts_per_seq[key], aa_variables_props_per_seq[key] = aa_variables_counts_and_props(aa_counts_per_seq[key])

        dict_for_files_aa_counts[file] = aa_counts_per_seq
        dict_for_files_aa_props[file] = aa_props_per_seq
        dict_for_files_variables_counts[file] = aa_variables_counts_per_seq
        dict_for_files_variables_props[file] = aa_variables_props_per_seq

    # One file per AA
    for aa in LAA:
        out = open(aa, 'w')
        out.write('Group,'+flsp[0:-1]+'\n')
        for group in dict_for_files_aa_counts.keys():
            count_of_aa = ''
            for specs in sorted(dict_for_files_aa_counts[group].keys()):
                count_of_aa += str(dict_for_files_aa_counts[group][specs][aa])+','
            out.write(group+',' + count_of_aa[0:-1] + '\n')
        out.close()
        os.system('mv %s %s/' %(aa, path_outputs_1))

    # one file per aa_variable
    for var in LV:
        out = open(var, 'w')
        out.write('Group,'+flsp[0:-1]+'\n')
        for group in dict_for_files_variables_counts.keys():
            count_of_var = ''
            for specs in sorted(dict_for_files_variables_counts[group].keys()):
                count_of_var += str(dict_for_files_variables_counts[group][specs][var])+','
            out.write(group+',' + count_of_var[0:-1] + '\n')
        out.close()
        os.system('mv %s %s/' %(var, path_outputs_2))

if __name__ == '__main__':
    main()