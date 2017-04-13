#!/usr/bin/python

## Author: Eric FONTANILLAS
## Date: 03.01.12
## Object: Recursively run the script "16_PATRON_stats.r" which compute stats and report graph pdf) for each branch reported in the codeml output


import string, os

## 1 ## Open input/output
branch_folder = "12_2_Branches"
L1 = os.listdir(branch_folder)

file_PATRON = open("16_PATRON_stats.r", "r")
PATRON = file_PATRON.read()

out_graph = "18_statsOutput_perBranches_GRAPHS"
out_table = "18_statsOutput_perBranches_TABLES"

for file in L1:
    print "\n\tProcess branch '%s' ..." %file[:-4]
    XXX_IN = "%s/%s" %(branch_folder, file)
    
    name_fasta = file[:-4]
    
    XXX_OUT_GRAPH = "%s/%s.pdf" %(out_graph, name_fasta)
    XXX_OUT_TABLE = "%s/%s.csv" %(out_table, name_fasta)

    
    New_PATRON = string.replace(PATRON, "XXX_IN", XXX_IN)
    New_PATRON = string.replace(New_PATRON, "XXX_OUT_GRAPH", XXX_OUT_GRAPH)
    New_PATRON = string.replace(New_PATRON, "XXX_OUT_TABLE", XXX_OUT_TABLE)

    print New_PATRON
    
    #print "\tOpen Patron ..."    
    PATRON_TMP = open("PATRON_tmp.r", "w")
    PATRON_TMP.write(New_PATRON)
    PATRON_TMP.close()
    
    os.system("R CMD BATCH PATRON_tmp.r")
    
    
