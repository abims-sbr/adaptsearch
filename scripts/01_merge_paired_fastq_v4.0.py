#!/usr/bin/python

## AUTHOR: Eric Fontanillas

## LAST VERSION: 14.04.2011

## DESCRIPTION: merge read1.fastq and read2.fastq files (from illumina 1.5 outputs with different numbers of read1 and read2)and remove sequences with low quality Phred score (if too much B at the end of a reads) 

import string, os, sys

#N=2    ## Fasta format
N=4    ## Fastq format

## CHECK IF ALL ARGUMENTS
if len(sys.argv) == 5:
    IN1 = open(sys.argv[1], "r")
    IN2 = open(sys.argv[2], "r")
    PREFIX = sys.argv[3]

else:
    print "\n\tThe command takes 4 arguments!! USAGE:\n\n\t\t$ merge_paired_fastq.py read1 read2 prefix_name_for_output Max_B_allowed\n"
    sys.exit()


## Phred quality score cutoff
B_cutoff = int(sys.argv[4])   ## Max of "B" authorized => if B_Cutoff > read length ===> means no cutoff are apply (no low quality sequence excluded!)



OUT1 = open("../tmp/%s_Bcutoff%d_paired.fas" %(PREFIX, B_cutoff), "w") ## paired fasta
OUT1_1 = open("../tmp/%s_Bcutoff%d_paired_read1.fas" %(PREFIX, B_cutoff), "w") ## paired fasta
OUT1_2 = open("../tmp/%s_Bcutoff%d_paired_read2.fas" %(PREFIX, B_cutoff), "w") ## paired fasta

OUT2 = open("../tmp/%s_Bcutoff%d_paired.fasq" %(PREFIX, B_cutoff), "w") ## paired fastq
OUT2_1 = open("../tmp/%s_Bcutoff%d_paired_read1.fasq" %(PREFIX, B_cutoff), "w") ## paired fastq
OUT2_2 = open("../tmp/%s_Bcutoff%d_paired_read2.fasq" %(PREFIX, B_cutoff), "w") ## paired fastq

OUT3 = open("../tmp/%s_Bcutoff%d_single1.fas" %(PREFIX, B_cutoff), "w") ## single 1 fasta
OUT4 = open("../tmp/%s_Bcutoff%d_single1.fasq" %(PREFIX, B_cutoff), "w") ## single 1 fastq
OUT5 = open("../tmp/%s_Bcutoff%d_single2.fas" %(PREFIX, B_cutoff), "w") ## single 2 fasta
OUT6 = open("../tmp/%s_Bcutoff%d_single2.fasq" %(PREFIX, B_cutoff), "w") ## single 2 fastq

OUT_read1_badQuality = open("../tmp/%s_Bcutoff%d_read1_badQuality.fasq" %(PREFIX, B_cutoff), "w")
OUT_read2_badQuality = open("../tmp/%s_Bcutoff%d_read2_badQuality.fasq" %(PREFIX, B_cutoff), "w")

LOG = open("../tmp/%s_Bcutoff%d_Counts_reads_and_pairs.log" %(PREFIX, B_cutoff), "w")

#i=1
j=1
b=0
bash1={}
while 1:
    new = IN1.readline()

    if not new:
        LOG.write("[Reads 1]\n")
        LOG.write("TOTAL reads number = %d reads\n" %j)
        percent_removed = float(b)/float(j)*100
        remaining = j - b
        percent_remaining = float(remaining)/float(j)*100
        LOG.write("Reads keeped (Good Phred Quality Score) = %d reads (%d %%)\n" %(remaining,percent_remaining))
        LOG.write("Reads removed (Too low Phred Quality Score) = %d reads (%d %%)\n\n" %(b,percent_removed))
        break

    if j%10000 ==0:
        print "[READ1] %d" %j

    if new[0] == "@":
        j=j+1
        name = new[1:-2]    ## remove "@" + "1" and "\n"
        name_fasta = ">" + new[1:-1]
        name_fastq = new[:-1]
        
        new = IN1.readline()
        seq = new[:-1]

        new = IN1.readline()  ## remove the name beginning with "+"
        phred_name = new[:-1]
        
        new = IN1.readline()  
        phred_quality = new[:-1]
        B = string.count(phred_quality, "B")
        
        if B <= B_cutoff:
            bash1[name] =[seq,name_fasta,name_fastq, phred_name, phred_quality]
        else:
            b=b+1
            OUT_read1_badQuality.write("%s\n" %name_fastq)
            OUT_read1_badQuality.write("%s\n" %seq)
            OUT_read1_badQuality.write("%s\n" %phred_name)
            OUT_read1_badQuality.write("%s\n" %phred_quality)

    #i=i+1

IN1.close()

#i=1
j=1
b=0
bash2={}
while 1:
    new = IN2.readline()

    if not new:
        #print "\n\tTOTAL Lines = %d lines\n" %i
        #seq_nb = j
        LOG.write("[Reads 2]\n")
        LOG.write("TOTAL reads number = %d reads\n" %j)
        percent_removed = float(b)/float(j)*100
        remaining = j - b
        percent_remaining = float(remaining)/float(j)*100
        LOG.write("Reads keeped (Good Phred Quality Score) = %d reads (%d %%)\n" %(remaining,percent_remaining))
        LOG.write("Reads removed (Low Phred Quality Score) = %d reads (%d %%)\n\n" %(b,percent_removed))
        break

    if j%10000 ==0:
        print "[READ2] %d" %j

    if new[0] == "@":
        j=j+1
        name = new[1:-2]    ## remove "@" + "1" and "\n"
        name_fasta = ">" + new[1:-1]
        name_fastq = new[:-1]
        
        new = IN2.readline()
        seq = new[:-1]

        new = IN2.readline()  ## remove the name beginning with "+"
        phred_name = new[:-1]
        
        new = IN2.readline()  
        phred_quality = new[:-1]
        B = string.count(phred_quality, "B")
        
        if B <= B_cutoff:
            bash2[name] =[seq,name_fasta,name_fastq, phred_name, phred_quality]
        else:
            b=b+1
            OUT_read2_badQuality.write("%s\n" %name_fastq)
            OUT_read2_badQuality.write("%s\n" %seq)
            OUT_read2_badQuality.write("%s\n" %phred_name)
            OUT_read2_badQuality.write("%s\n" %phred_quality)

    #i=i+1

IN2.close()

######################################
p = 0
s1 = 0
s2 = 0
for name in bash1.keys():
    list1 = bash1[name]
    seq1 = list1[0]
    name_fasta1 = list1[1]
    name_fastq1 = list1[2]
    name_phred1 = list1[3]
    phred_quality1 = list1[4]

    ## 1 ## Save "paired 1 and 2 reads"
    try:
        list2 = bash2[name]
        seq2 = list2[0]
        name_fasta2 = list2[1]
        name_fastq2 = list2[2]
        name_phred2 = list2[3]
        phred_quality2 = list2[4]
    
        ## a ## Print fasta output
        OUT1.write("%s\n" %name_fasta1)
        OUT1.write("%s\n" %seq1)
        OUT1_1.write("%s\n" %name_fasta1)
        OUT1_1.write("%s\n" %seq1)

        OUT1.write("%s\n" %name_fasta2)
        OUT1.write("%s\n" %seq2)
        OUT1_2.write("%s\n" %name_fasta2)
        OUT1_2.write("%s\n" %seq2)

        ## b ## Print fastq output
        p=p+1
        OUT2.write("%s\n" %name_fastq1)
        OUT2.write("%s\n" %seq1)
        OUT2.write("%s\n" %name_phred1)
        OUT2.write("%s\n" %phred_quality1)
        OUT2_1.write("%s\n" %name_fastq1)
        OUT2_1.write("%s\n" %seq1)
        OUT2_1.write("%s\n" %name_phred1)
        OUT2_1.write("%s\n" %phred_quality1)

        
        OUT2.write("%s\n" %name_fastq2)
        OUT2.write("%s\n" %seq2)
        OUT2.write("%s\n" %name_phred2)
        OUT2.write("%s\n" %phred_quality2)

        OUT2_2.write("%s\n" %name_fastq2)
        OUT2_2.write("%s\n" %seq2)
        OUT2_2.write("%s\n" %name_phred2)
        OUT2_2.write("%s\n" %phred_quality2)
                        

        
    ## 2 ## Save "single 1 reads"
    except:
        ## a ## Print fasta output
        OUT3.write("%s\n" %name_fasta1)
        OUT3.write("%s\n" %seq1)

        ## b ## Print fastq output
        s1 = s1 + 1
        OUT4.write("%s\n" %name_fasta1)
        OUT4.write("%s\n" %seq1)
        OUT4.write("%s\n" %name_phred1)
        OUT4.write("%s\n" %phred_quality1)
        
        continue
        
    

for name in bash2.keys():
    list2 = bash2[name]
    seq2 = list2[0]
    name_fasta2 = list2[1]
    name_fastq2 = list2[2]
    name_phred2 = list2[3]
    phred_quality2 = list2[4]

    try:  ## Already done
        list1 = bash1[name]
        seq1 = list1[0]
        name_fasta1 = list1[1]
        name_fastq1 = list1[2]
        name_phred1 = list1[3]
        phred_quality1 = list1[4]

    except: ## 3 ## Save "single 2 reads"
        ## a ## Print fasta output
        OUT5.write("%s\n" %name_fasta2)
        OUT5.write("%s\n" %seq2)

        ## b ## Print fastq output
        s2 = s2 + 2
        OUT6.write("%s\n" %name_fasta2)
        OUT6.write("%s\n" %seq2)
        OUT6.write("%s\n" %name_phred2)
        OUT6.write("%s\n" %phred_quality2)

        continue

LOG.write("[Merged]\n")
LOG.write("Pairs recovered = %d\n" %p)
LOG.write("Short single from read1 = %d\n" %s1)
LOG.write("Short single from read2 = %d\n" %s2)



OUT1.close()
OUT2.close()
OUT3.close()
OUT4.close()
OUT5.close()
OUT6.close()
OUT_read1_badQuality.close()
OUT_read2_badQuality.close()
