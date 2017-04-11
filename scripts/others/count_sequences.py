#!/usr/bin/python

## AUTHOR: Eric Fontanillas

## LAST VERSION: 14.04.2011

## DESCRIPTION: count lines

N=2    ## Fasta format
#N=4    ## Fastq format

import string, os, sys

path_IN = sys.argv[1]

IN = open(path_IN, "r")


i=1
while 1:
    new = IN.readline()

    if i%100000 ==0:
        print i
    
    if not new:
        print "\n\tTOTAL Lines = %d lines\n" %i
        seq_nb = (i-1)/N
        print "\tTOTAL sequence number = %d reads\n" %seq_nb
        break
    
    i=i+1
