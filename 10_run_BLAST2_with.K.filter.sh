#!/bin/bash

INPUT_DB="$1"
INPUT_QUERY="./09_onlyMatch_filtered_nucleotideBACK.fasta"

OUTPUT="./11_outputBlast.txt"
LOG="./11_log"

## 1.1. Format database ##
# PARAMETERS:
# -i : Database to format
# -p : input type : T: proteins // F: nucleotides
# -o : Parses deflines and indexes seqIDs

echo '***** START format Database [BLAST] *****'
formatdb -i  $INPUT_DB -p F -o T
echo '***** END format Database [BLAST] *****'

### 1.2. Run Blast ###
# PARAMETERS:
# -p : program: "blastp", "blastn", "blastx", "tblastn", or "tblastx" 
# -d : database
# -i : query
# -e : Expectation value => -e 1e-10 [default = 0]
# -F : DUST filter => -F T/F [default = T]
# -T : Produce HTML output => -T T/F [default = F]

echo '***** START run BLAST *****'
time blastall -p tblastx -d $INPUT_DB -i $INPUT_QUERY -o $OUTPUT -T F -e  1e-20 -F "m S" -b 1 -v 1 -K 1 > $LOG
echo '***** END run BLAST *****'
echo ' '

####################
### 2. Parse Blast Output ###
####################
# Parse Blast Output
# Produce output with Pairwise Matches
# INPUT = "05_blast_ApoSouth_vs_ApoNorth.txt"
# OUTPUT = "07_PairwiseMatch.fasta"
# LOG = "08_parseBLASToutput.log"

#echo '***** START parsing of the Blast Output *****'
#./03_scriptExtractMatch_v9.py
#echo '***** END parsing of the Blast Output *****'




##############################



