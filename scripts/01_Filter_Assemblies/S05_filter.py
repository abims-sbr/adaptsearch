#!/usr/bin/env python
#filters the sequences depending on their length after cap3, makes the sequences names compatible with the phylogeny workflow
#python filter.py file length_threshold_nucleotides output

import string, os, sys, re

path_IN = sys.argv[1]
threshold=int(sys.argv[2]) #minimum number of nucleotides for one sequence
file_OUT = open(sys.argv[3], "w")
f_in = open(path_IN, "r")
inc=1
while 1:
    line=f_in.readline()
    if not line:
        break
    line=f_in.readline()
    name=">"+path_IN[:2]+str(inc)+"_1/1_1.000_"
    if len(line)-1>threshold-1:
        inc+=1
        file_OUT.write("%s" %name)
        file_OUT.write(str(len(line)-1)+"\n")
        file_OUT.write("%s" %line)
f_in.close()
file_OUT.close()
    
    

#filtre eventuel sur les petits transcrits





  
