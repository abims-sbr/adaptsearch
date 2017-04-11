#!/usr/bin/python

# formatting a fasta format into phylip format for using with PAML

import string, os, sys


path_IN = sys.argv[1]
path_OUT = sys.argv[2]


###### def1 ######
def extract_interleavedFASTA(path_IN, path_OUT):
      file_IN = open(path_IN, "r")
      file_OUT = open(path_OUT, "w")

      i = 0
      n = 1
      while 1:
            nextline = file_IN.readline()
            if not nextline:
                  #####################################################
                  ### III. ### LAST SEQUENCE [Only record sequence] ###
                  #####################################################
                  if i != 0:
                        file_OUT.write("%s\n" %fasta_name[:-1])
                        new_fasta_seq = string.replace(fasta_seq, "\n", "")
                        file_OUT.write("%s\n" %new_fasta_seq)
                        print n
                  break

            if nextline[0] == ">":
                  ########################################################
                  ### I. ### FIRST SEQUENCE [Only initialize sequence] ###
                  ########################################################
                  if i == 0:                   
                        fasta_name = nextline
                        fasta_seq = "" # initialise the fasta sequence
                        i = 1  ## Will indicates that the first sequence is treated
                  #####################################################################################
                  ### II. ### OTHER SEQUENCES [Record previous sequence + Initialize next sequence] ###
                  #####################################################################################
                  else:                            
                        ## 1 ## Record the previous sequence:
                        file_OUT.write("%s\n" %fasta_name[:-1])
                        new_fasta_seq = string.replace(fasta_seq, "\n", "")
                        file_OUT.write("%s\n" %new_fasta_seq)
                        print n

                        ## 2 ## Initialize the next sequence:
                        fasta_name = nextline
                        fasta_seq = ""

                        n = n+1
                        
            else:
                  fasta_seq = fasta_seq + nextline
            

      #dico = {}
      #      file_txt = ""
      
      #       subfile = file.read()
            
      #       L1 = string.split(subfile, '>')

      #       for element in L1:

      #           if element != '':
      #               j = string.find(element, '\n')
      
      #               fasta_name = element[:j]
      
      
      #               seq = element[j+1:-1]
      #               sequence = string.replace(seq, '\n', '')
      
      #               file_txt = file_txt + ">%s\n" %fasta_name + "%s\n"   %sequence
      file_IN.close()
      file_OUT.close()
#-#-#-#-#-#-#-#-#-#-#

###################
### RUN RUN RUN ###
###################
import string
file = extract_interleavedFASTA(path_IN, path_OUT)   ### DEF1 ###
