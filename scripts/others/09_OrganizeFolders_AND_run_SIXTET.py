#!/usr/bin/python

## Author: Eric FONTANILLAS
## Date: 06.12.12
## Object: Organize folder for run with the 6 species (Ap, Ac, Pp, Pg, Ps, Pf)


##################
###### DEF1 ######
##################
# Dans un multialignement fasta, cette fonction permet de formatter les noms de chaque sequence fasta dans le cas de COMPARAISON DE SEQUENCE PAR PAIR!!

def Get_SIXTET(F, number_sp):
    list_locus = [] 
    lists = []

    list_SIXTET = F.readlines()
    list_SIXTET = "".join(list_SIXTET)
    list_SIXTET = list_SIXTET.split("\n")
    del list_SIXTET[-1]

    number = number_sp*2
    i=0
    j=0
    m=0
    z=0
    while z < len(list_SIXTET) :
	if list_SIXTET[z][0] == ">" :
	    list_SIXTET[z] = list_SIXTET[z].replace(list_SIXTET[z][0], "")
	z+=1	

    while m < len(list_SIXTET)+2 :
	i+=1
        if i <= number :
            list_locus.append(list_SIXTET[j])
        else:
            lists.append(list_locus)
	    list_locus = []
	    j = j-1
	    i=0
	m+=1
        j+=1

    return(lists)
##############################################


###################
### RUN RUN RUN ###
###################

import string, os, zipfile, sys, re

## 1 ## Open input/output

## create a file with all loci
os.mkdir("./locus/")
zfile= zipfile.ZipFile(sys.argv[1])
for name in zfile.namelist() :
    zfile.extract(name, "./locus/")
path_IN1 = "./locus/"
L_IN = os.listdir(path_IN1)

F1 = open("08_All_Loci_in_1_flat_file.fasta", "w")
for files in L_IN :
    F0 = open("%s/%s" %(path_IN1,files), "r")
    next = F0.readlines()
    number_seq = 0
    for seq in next :
	number_seq += 1
	F1.write(seq)
    F0.close()
F1.close
number_sp = number_seq / 2 # number of species
os.system("rm -rf locus/")
F1 = open("08_All_Loci_in_1_flat_file.fasta", "r")

# create a file with a patron tree with the output of RAxML
i = 0
numbers = "[0-9]"
F3 = open("00_PATRON_TREE_FreeRatios.tree", "w")
F0 = open(sys.argv[2], "r")
next = F0.read()
next = next.split(":")
next = "".join(next)
while i<len(next) :
    if re.match(numbers, next[i]) or re.match("\.", next[i]) :
	next = next.replace(next[i], "")
    i+=1
print next
F3.write(next)
F0.close
F3.close

F3 = open("00_PATRON_TREE_FreeRatios.tree", "r")
tree_patron_null = F3.read()
F3.close()
F3 = open("00_PATRON_TREE_FreeRatios.tree", "r")
tree_patron_alternate = F3.read()
F3.close()

os.mkdir("./TREE_FreeRatios")
path = "./TREE_FreeRatios"

## Read fasta file and split the multilocus file in several files
list_SIXTET = Get_SIXTET(F1, number_sp)    ### DEF 1 ###

d = 1
for sixtet in list_SIXTET:
    file_name = "locus_%d" %d
    d = d+1
    subpath_null = "%s/%s_NULL" %(path, file_name)
    os.mkdir(subpath_null)
    subpath_alternate = "%s/%s_ALTERNATE" %(path, file_name)
    os.mkdir(subpath_alternate)
       
    ## copy codeml (from patron)
    # NULL
    os.system("cp 00_PATRON_codeml_NULL.ctl %s/codeml.ctl" %subpath_null)
    # ALTERNATE
    os.system("cp 00_PATRON_codeml_ALTERNATE.ctl %s/codeml.ctl" %subpath_alternate)

    sixtet_fasta_null = open("%s/sixtet.fasta" %subpath_null, "w")
    sixtet_fasta_alternate = open("%s/sixtet.fasta" %subpath_alternate, "w")
    new_tree_patron_null = tree_patron_null
    tree_file_null = open("%s/sixtet.tree" %subpath_null, "w")
    new_tree_patron_alternate = tree_patron_alternate
    tree_file_alternate = open("%s/sixtet.tree" %subpath_alternate, "w")

    k = 0
    n = k+1
    while n<len(sixtet)+1 :
	name = sixtet[k]
	print name
	short_name = name[0:2]
	seq = sixtet[n]

	ex1 = "(%s" %short_name
	new_ex1 = "(%s" %name
	ex2 = "%s)" %short_name
	new_ex2 = "%s)" %name

        new_tree_patron_null = new_tree_patron_null.replace(ex1, new_ex1, 1)
	new_tree_patron_null = new_tree_patron_null.replace(ex2, new_ex2, 1)
	new_tree_patron_alternate = new_tree_patron_alternate.replace(ex1, new_ex1, 1)
	new_tree_patron_alternate = new_tree_patron_alternate.replace(ex2, new_ex2, 1)

	sixtet_fasta_null.write(">%s\n" %name)
	sixtet_fasta_null.write("%s\n" %seq)
	sixtet_fasta_alternate.write(">%s\n" %name)
	sixtet_fasta_alternate.write("%s\n" %seq)

	k+=2
	n+=2

    sixtet_fasta_null.close()
    sixtet_fasta_alternate.close()

    tree_file_null.write(new_tree_patron_null)
    tree_file_null.close()
    tree_file_alternate.write(new_tree_patron_alternate)
    tree_file_alternate.close()

    ## convert alignment file from fasta to phylip AND delete input fasta file
    # NULL
    CURRENT_PATH = os.path.abspath("./")  ## Allow to go back to the script path
    
    os.system("cp 00_PATRON_fasta2phylip.py %s/fasta2phylip.py" %subpath_null)
    os.chdir("%s" %subpath_null)
    os.system("./fasta2phylip.py")
    os.system("codeml codeml.ctl")  ## execute codeml
    os.chdir(CURRENT_PATH)

    # ALTERNATE
    os.system("cp 00_PATRON_fasta2phylip.py %s/fasta2phylip.py" %subpath_alternate)
    os.chdir("%s" %subpath_alternate)
    os.system("./fasta2phylip.py")
    os.system("codeml codeml.ctl")  ## execute codeml
    os.chdir(CURRENT_PATH)
    
F1.close()
F3.close()
