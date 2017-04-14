#!/usr/bin/env python
import sys, string, os, itertools, re, zipfile

# coding: utf8

i=1
j=i+1
L1 = []
fas = "^.*fas$"
fasta = "^.*fasta$"
script_path = os.path.dirname(sys.argv[0])

infiles = sys.argv[1]
assembler = sys.argv[2]
length_seq_max = sys.argv[3]
percent_identity = sys.argv[4]
overlap_length = sys.argv[5]

for name in str.split(infiles,","):
    suffix=name[:2]

    #Fasta sequences on one line
    name_fasta_formatter = "01%s" %name
    os.system("cat '%s' | fasta_formatter -w 0 -o '%s'" % (name,name_fasta_formatter))

    if assembler=="velvet" :
        name_remove_redondancy="02_%s" %name
        os.system("python %s/S02a_remove_redondancy_from_velvet_oases.py %s %s" %(script_path,name_fasta_formatter, name_remove_redondancy))

        name_find_orf_input="%s%s" %(suffix,name)
        os.system("sed -e 's/Locus_/%s/g' -e 's/_Confidence_/_/g' -e 's/_Transcript_/_/g' -e 's/_Length_/_/g' %s > %s" % (suffix,name_remove_redondancy,name_find_orf_input))

    elif assembler=="trinity" :
        #replace white space by "_"
        name_sed="02%s" %name
        os.system("sed -e 's/ /_/g' %s > %s " % (name_fasta_formatter,name_sed))

        #Format the name of the sequences with good name
        name_format_fasta="03%s" %name
        os.system("python %s/S02b_format_fasta_name_trinity.py %s %s %s" %(script_path,name_fasta_formatter, name_format_fasta, suffix))

        #Apply first script to avoid reductant sequences
        name_find_orf_input="04%s" %name
        os.system("python %s/S03_choose_one_variants_per_locus_trinity.py %s %s" %(script_path,name_format_fasta, name_find_orf_input))

    #Pierre guillaume find_orf script for keeping the longuest ORF
    name_find_orf="05%s"% name
    os.system("python %s/S04_find_orf.py %s %s" %(script_path,name_find_orf_input,name_find_orf))
    #Apply cap3
    os.system("cap3  %s -p %s -o %s"%(name_find_orf,percent_identity,overlap_length))
    #Fasta sequences on one line
    #Il faudrait faire un merge des singlets et contigs! TODO
    os.system("zcat -f < '%s.cap.singlets' | fasta_formatter -w 0 -o '%s'" % (name_find_orf,suffix))
    #Apply pgbrun script filter script TODO length parameter
    name_filter="%s%s"%(suffix,name)
    os.system("python %s/S05_filter.py %s %s %s" %(script_path,suffix,length_seq_max,name_filter))
    L1.append(name_filter)
    f = zipfile.ZipFile("sequences_filtered.zip", "w")
    for files in L1 :
        f.write("./%s" %files)
