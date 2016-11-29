#!/usr/bin/python

## AUTHOR: Eric Fontanillas
## LAST VERSION: 20/08/14 by Julie BAFFARD

## DESCRIPTION: Prepare to run multialign on assemblages on several cluster nodes

import os,sys
script_path = os.path.dirname(sys.argv[0])
###############################################
### DEF 1 : Split a list in several sublist ###
###############################################
def chunks(list, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(list), n):
        yield list[i:i+n]
######################################


##########################################
### DEF 2 : Prepare run for blastalign ###
##########################################
def prepare_BLASTALIGN_runs(list_file):

    ln = len(list_file)
    i = 0
    list_of_sublist = list(chunks(list_file, 5000))           ### DEF2 ###

    list_files_failed = []

    k=0
    for sublist in list_of_sublist:
        for fasta_file in sublist:
            i = i + 1
            S1 = string.split(fasta_file, ".")
            fasta_name = S1[0]

	    # filter the "N"
	    f = open("./%s.fasta" %fasta_name, "r")
	    nextline=f.readlines()
	    f.close()

	    j = 0
	    while j<len(nextline) :
		if not nextline[j].startswith(">") :
		    nextline[j] = nextline[j].upper()
		    nombre = nextline[j].rfind("N")
		    if nombre != -1 :
			nextline[j] = nextline[j].replace("N", "")
		j+=1
	    nextline = "".join(nextline)

	    files = open("./%s.fasta" %fasta_name, "w")
	    files.write(nextline)
	    files.close()

            ## run individual script
	    os.system("%s/BlastAlign -m %s -n %s -i ./%s.fasta\n" %(script_path,sys.argv[3],sys.argv[4],fasta_name))

	    try:
        	phylip_file = open("./%s.fasta.phy" %fasta_name, "r")       
    	    except IOError:
		list_files_failed.append(fasta_file)

	    if sys.argv[2] == "oui" :
                os.system("python %s/phylip2fasta.py ./%s.fasta.phy ./%s.fasta.fasta\n" %(script_path,fasta_name, fasta_name))
		os.system("rm -f ./%s.fasta\n" %fasta_name)
		os.system("mv ./%s.fasta.fasta ./%s.fasta\n" %(fasta_name, fasta_name))
  
    return(list_files_failed)
            
######################################

###################
### RUN RUN RUN ###
###################

import string, os, time, re, sys, zipfile, re

## 1 ## INPUT/OUTPUT
list_file = []
zfile = zipfile.ZipFile(sys.argv[1])
for name in zfile.namelist() :
    list_file.append(name)
    zfile.extract(name, "./")

## 2 ## RUN
list_files_failed = prepare_BLASTALIGN_runs(list_file)   ### DEF2 ###
f_failed = open("./list_files_failed.txt", "w")
f_failed.write("Number of files failed with BlastAlign : %d\n" %len(list_files_failed))
for files in list_files_failed :
    f_failed.write("\t%s \n" %files)


#Convertion in zip format
f_phylip = zipfile.ZipFile("Alignment_locus_phy.zip", "w")
f_nexus = zipfile.ZipFile("Alignment_locus_nxs.zip", "w")
f_fasta = zipfile.ZipFile("Alignment_locus_fasta.zip", "w")

phylip = "^.*fasta.phy$"
nexus = "^.*fasta.nxs$"
fasta = "^.*fasta$"

folder = os.listdir("./")

for files in folder :
    if re.match(phylip, files) :
	f_phylip.write("./%s" %files)
    if re.match(nexus, files) :
	f_nexus.write("./%s" %files)
    if re.match(fasta, files) :
	f_fasta.write("./%s" %files)
