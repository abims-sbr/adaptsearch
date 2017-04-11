#!/usr/bin/python

# formatting a fasta format into phylip format for using with PAML

import string, os

# while 1:
#     f1 = raw_input("\nFichier FASTA a ouvrir : ")

#     try:
#         F1 = open("%s" % f1, 'r')
#         break
#     except IOError :
#         print "\nLe fichier n'existe pas !!"
#         continue

F1 = open("sixtet.fasta", "r")

F2 = open("sixtet.phy", 'w')

###### def1 ######
# Dans un multialignement fasta, cette fonction permet de formatter les noms de chaque sequence fasta

def format(F):
    c = 0
    fichier = ""
    while 1 :
        next = F.readline()
        if not next :
            break
        if next[0] == ">" :
            if c > 0 : # add "\n" before each ">" except for the 1rst line
                fichier = fichier + "\n"
            #fichier = fichier + next[1:-1] + "\n"
            fichier = fichier + next[1:]

            c = c + 1

        else :
            newline = next[:-1]
            fichier = fichier + newline
    return (fichier,c)
#-#-#-#-#-#-#-#-#-#-#


list = format(F1)  ## def1 ##
newfile = list[0]
taxaNumber = list[1]

F3 = open("temp_file" , 'w')
F3.write(newfile)
F3.close()

F4 = open("temp_file", 'r')


next = F4.readline() # read 1rst line
next = F4.readline() # read 2nd line
seq_length = len(next[:-1]) # [:-1] car sinon je compte "\n" en plus dans la longueur de la sequence !!!

F2.write("   %d %d\n" % (taxaNumber, seq_length))
F2.write(newfile)

F4.close()
os.system("rm temp_file")

F1.close()
F2.close()
