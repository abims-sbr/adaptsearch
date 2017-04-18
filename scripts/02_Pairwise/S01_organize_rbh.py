#!/usr/bin/env python
## AUTHOR: Eric Fontanillas
## LAST VERSION: 03/09/14 by Julie BAFFARD

import sys, string, os, itertools, re, zipfile
from multiprocessing import Pool

infiles = sys.argv[1]
evalue = sys.argv[2]
thread = sys.argv[3]

def launching(RBH_folder):
    os.system("bash %s/XXX_pipeline_%s.sh %s" %(RBH_folder,RBH_folder,RBH_folder))

L2 = str.split(infiles,",") #list on input file
list_RBH = []
m=1

list_pairwise = itertools.combinations(L2, 2) # Sets couple sp1/sp2

for pairwise in list_pairwise :
    print "Pair of species:"
    print pairwise
    # simplifies the file names (ex : remove "_Trinity.fasta")
    DIR1 = pairwise[0]
    DIR2 = pairwise[1]

    S1 = string.split(DIR1, ".")
    S1 = S1[0]
    S1 = string.split(S1, '_')
    shortDIR1 = S1[0]

    S2 = string.split(DIR2, ".")
    S2 = S2[0]
    S2 = string.split(S2, "_")
    shortDIR2 = S2[0]

    RBH_folder = shortDIR1 + "_" + shortDIR2
    list_RBH.append(RBH_folder) # List of all couples with simplified names

    os.mkdir("./%s" %RBH_folder)

    os.system("cp -fr S03_run_blast_with_k_filter.sh %s/" %(RBH_folder))
    os.system("cp -fr S04_run_blast2_with_k_filter.sh %s/" %(RBH_folder))

    if L2 != [] :
        os.system("cp -fr %s ./%s/%s" %(DIR1, RBH_folder, DIR1))
    os.system("cp -fr %s ./%s/%s" %(DIR2, RBH_folder, DIR2))

    pipeline_script = open("./%s/XXX_pipeline_%s.sh" %(RBH_folder,RBH_folder), "w")
    #pipeline_patron = open("%s/XXX_patronPipeline.sh" %(SCRIPTPATH), "r")
    pipeline_patron = open("S02_xxx_patron_pipeline.sh", "r")

    #pipeline_patron = open("/home/umr7144/abice/vmataigne/Documents/AdaptSearch/adaptsearch-master/scripts/XXX_patronPipeline.sh", "r")

    # Setting parameters in XXX_Patron_Pipeline
    F1 = pipeline_patron.read()
    F1 = string.replace(F1, "XX1", "./%s/%s" %(RBH_folder,DIR1))
    F1 = string.replace(F1, "XX2", "./%s/%s" %(RBH_folder,DIR2))
    F1 = string.replace(F1, "XX3", DIR1)
    F1 = string.replace(F1, "XX4", DIR2)
    F1 = string.replace(F1, "XX5", evalue)
    F1 = string.replace(F1, "WORK_DIR", "%s" %RBH_folder)

    pipeline_script.write(F1)

    pipeline_patron.close()
    pipeline_script.close()

## Multithreading
pool = Pool(processes=int(thread))
result = pool.map(launching, list_RBH)

os.system("rm 25_DNA*")
os.system("rm 19_Reci*")
os.system("rm *.fasta")
