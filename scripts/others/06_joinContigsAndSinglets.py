#!/usr/bin/python

# formatting a fasta format into phylip format for using with PAML

import string, os, sys


path_IN1 = sys.argv[1]    # contigs
path_IN2 = sys.argv[2]    # singlets
path_OUT = sys.argv[3]


IN1 = open(path_IN1, "r")
IN2 = open(path_IN2, "r")
OUT = open(path_OUT, "w")

while 1:
      next = IN1.readline()
      if not next:
            break
      next = string.replace(next, "Contig", "ContigI_")   ## "I" : cause 1rst run of CAP3
      OUT.write(next)
while 1:
      nexteu = IN2.readline()
      if not nexteu:
            break
      OUT.write(nexteu)

IN1.close()
IN2.close()
OUT.close()

