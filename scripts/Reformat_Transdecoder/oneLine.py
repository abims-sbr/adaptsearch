#!/usr/bin/env python

import sys, string
from Bio import SeqIO # BioPython

def seqOneLine(file, name): #file : $2 (temp file) in script.sh
    n = file.split('.')    
    #name = "%s_oneline.fasta" %n[0]   
    
    new_file = open(name, "w")
    with new_file:
        for seq_record in SeqIO.parse(file, "fasta"):
            gid = seq_record.id
            gid = gid.split("(")
            gid = gid[0]
            new_file.write(gid)
            new_file.write("\n")
            new_file.write(str(seq_record.seq))
            new_file.write("\n")

def main():
    seqOneLine(sys.argv[1], sys.argv[2])

if __name__ == "__main__":
    main()
