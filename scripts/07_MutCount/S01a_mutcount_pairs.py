#!/usr/bin/env python
# coding: utf-8
#

import os, itertools,sys


list_species=[]
num_sampled= sys.argv[1]
num_iter=sys.argv[2]
list_species= sys.argv[3]
list_species=list_species.split(",")
pairs_list=list(itertools.combinations(list_species,2))
print pairs_list

pairs_file=open("pairs.txt","w")
pairs_file.write("background: %s %s + all\n"%(num_sampled,num_iter))
for pair in pairs_list:
	pairs_file.write("%s %s\n"%(pair[0],pair[1]))
pairs_file.close()





