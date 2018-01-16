#!/usr/bin/env python
#coding: utf-8

## AUTHOR: Eric Fontanillas
## LAST VERSION: 10.2017 by Victor Mataigne

import glob, sys, string, os
from Bio import SeqIO # BioPython

# Improvments to do:
    # code factoring if possible
    # simplify the method (too many intermediate files with weird names)
    
def seqOneLine(file):
    n = file.split('.')    
    name = "%s_oneline.fasta" %n[0]   
    
    new_file = open(name, "w")
    with new_file:
        for seq_record in SeqIO.parse(file, "fasta"):
            gid = seq_record.id
            gid = gid.split("(")
            gid = gid[0]
            new_file.write(">")
            new_file.write(gid)
            new_file.write("\n")
            new_file.write(str(seq_record.seq))
            new_file.write("\n")

def nameFormatting(name, script_path, prefix):
    f = open(name, "r")    
    f1 = f.readline() # Only need to check first line to know the assembler which has been used
    f.close()
    name_find_orf_input = ""

    if f1.startswith(">Locus"):
        name_remove_redondancy = "02_%s" %name
        os.system("python S02a_remove_redondancy_from_velvet_oases.py %s %s" %(name, name_remove_redondancy))
        name_find_orf_input = "%s%s" %(prefix, name)
        os.system("sed -e 's/Locus_/%s/g' -e 's/_Confidence_/_/g' -e 's/_Transcript_/_/g' -e 's/_Length_/_/g' %s > %s" % (prefix, name_remove_redondancy, name_find_orf_input))
    elif f1.startswith(">c"):        
        #Format the name of the sequences with good name
        name_format_fasta = "03%s" %name
        os.system("python S02b_format_fasta_name_trinity.py %s %s %s" %(name, name_format_fasta, prefix))
        #Apply first script to avoid reductant sequences
        name_find_orf_input = "04%s" %name
        os.system("python S03_choose_one_variants_per_locus_trinity.py %s %s" %(name_format_fasta, name_find_orf_input))

    return name_find_orf_input

def main():
    os.mkdir("outputs")
    script_path = os.path.dirname(sys.argv[0])    
    length_seq_max = sys.argv[2]
    percent_identity = sys.argv[3]
    overlap_length = sys.argv[4]

    #for name in str.split(sys.argv[1], ","):
        #seqOneLine(name)

    #path = glob.glob('*_oneline.fasta')    

    #for name in path:
    for name in str.split(sys.argv[1], ","):         
        prefix=name[0:2]

        #debut modif
        name_fasta_formatter = "01%s" %name
        os.system("cat '%s' | fasta_formatter -w 0 -o '%s'" % (name, name_fasta_formatter))
        name_find_orf_input = nameFormatting(name_fasta_formatter, script_path, prefix)
        # fin modif

        #name_find_orf_input = nameFormatting(name, script_path, prefix)
        #Pierre guillaume find_orf script for keeping the longuest ORF
        name_find_orf = "05%s"% name
        os.system("python S04_find_orf.py %s %s" %(name_find_orf_input, name_find_orf))
        #Apply cap3
        os.system("cap3 %s -p %s -o %s"%(name_find_orf, percent_identity, overlap_length))
        #Il faudrait faire un merge des singlets et contigs! TODO
        os.system("zcat -f < '%s.cap.singlets' | fasta_formatter -w 0 -o '%s'" % (name_find_orf, prefix))
        #Apply pgbrun script filter script TODO length parameter
        name_filter = "%s%s"%(prefix, name)
        os.system("python S05_filter.py %s %s outputs/%s" %(prefix, length_seq_max, name_filter))

        #path = glob.glob('outputs/*_oneline.fasta')
        
    #for name in path:
        #name_new = name.split("_oneline.")[0]
        #os.system("mv %s %s.fasta" %(name, name_new))

if __name__ == "__main__":
    main()
