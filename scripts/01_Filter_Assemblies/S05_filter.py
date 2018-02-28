#!/usr/bin/env python
#filters the sequences depending on their length after cap3, makes the sequences names compatible with the phylogeny workflow
#python filter.py file length_threshold_nucleotides output

import string, os, sys, re, itertools

path_IN = sys.argv[1]
threshold = int(sys.argv[2]) #minimum number of nucleotides for one sequence
file_OUT = open(sys.argv[3], "w")
inc = 1
with open(path_IN, "r") as f_in:
    for ignored, sequence in itertools.izip_longest(*[f_in]*2):
        name=">"+path_IN[:2]+str(inc)+"_1/1_1.000_"
        if len(sequence)-1>threshold-1:
            inc+=1
            file_OUT.write("%s" %name)
            file_OUT.write(str(len(sequence)-1)+"\n")
            file_OUT.write("%s" %sequence)
file_OUT.close()

#filtre eventuel sur les petits transcrits