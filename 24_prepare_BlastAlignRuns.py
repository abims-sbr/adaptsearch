#!/usr/bin/python

## AUTHOR: Eric Fontanillas

## LAST VERSION: 09.08.2010

## DESCRIPTION: Prepare to run multialign on assemblages on several cluster nodes


##################################
### DEF 1. Clean output folder ###
##################################
def clean_folder(path_folderOUT):
    print "Cleaning of output folder in progress ..."
    
    
    L_test_path_empty = os.listdir(path_folderOUT)

    if L_test_path_empty != []:
        os.system("rm -fr %s/*" %path_folderOUT)
    print "Cleaning done!"

#########################################
### DEF 2. Prepare run for blastalign ###
#########################################
def prepare_BLASTALIGN_runs(path_folderIN, path_folderOUT, path_BLASTALIGN_BINARIES):

    L_IN = os.listdir(path_folderIN)

    ln = len(L_IN)

    script_sh = open("26_general_script_runs_blastAlign.sh", "w")
    
    i = 0
    for fasta_file in L_IN:
        i = i + 1
        print "%d/%d" %(i,ln)
        S1 = string.split(fasta_file, "_")
        fasta_name = S1[0]

        path_fileIN = "%s/%s" %(path_folderIN, fasta_file)
        path_subFolderOUT = "%s/%s" %(path_folderOUT, fasta_name)

        path_fileOUT = "%s/%s.fasta" %(path_subFolderOUT, fasta_name)

        ## 1 ## create subfolder
        os.mkdir(path_subFolderOUT)

        ## 2 ## copy fasta file (unaligned trio) in the subfolder
        os.system("cp %s %s" %(path_fileIN,  path_fileOUT))

        ## 3 ## copy BLASTALIGN binaries in the subfolder
        os.system("cp %s/* %s/" %(path_BLASTALIGN_BINARIES, path_subFolderOUT))

        ## 4 ## write individual script

        script = open("%s/script_eric.sh" %path_subFolderOUT, "w")
        script.write("BlastAlign -i %s.fasta\n" %fasta_name)
        script.write("phylip2fasta.py %s.fasta.phy %s.fasta.fasta\n" %(fasta_name,fasta_name))
        os.system("chmod 755 %s/script_eric.sh" %path_subFolderOUT)
        script.close()

        ## 5 ## write general script
        script_sh.write("cd %s\n" %path_subFolderOUT)
        script_sh.write("./script_eric.sh\n")
        script_sh.write("cd ../../\n")

    script_sh.close()
    return()
            
######################################



######################################
### RUN RUN RUN ###
###################

import string, os, time, re, sys

path_folderIN = "21_TRIO_UNALIGNED"
path_folderOUT = "25_run_BLASTALIGN"
path_BLASTALIGN_BINARIES = "23_BlastAlign_binary"

## 1 ## Clean folder from previous runs
clean_folder(path_folderOUT)    ### DEF1 ###

## 2 ## Prepare run for blastalign
prepare_BLASTALIGN_runs(path_folderIN, path_folderOUT, path_BLASTALIGN_BINARIES)   ### DEF2 ###

