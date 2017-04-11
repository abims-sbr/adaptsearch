#!/bin/env bash

SHORT_FILE=$1 #short-name-query_short-name-db

INPUT_DB=$2
SHORT_DB=$4 

INPUT_QUERY="./$SHORT_FILE/09_onlyMatch_filtered_nucleotideBACK_$SHORT_FILE.fasta"
SHORT_QUERY=$3

option_e=$5 

OUTPUT="./$SHORT_FILE/11_outputBlast_$SHORT_FILE.txt"
LOG="./$SHORT_FILE/11_$SHORT_FILE.log" 

## 1.1. Format database ##
# PARAMETERS:
# -i : Database to format
# -p : input type : T: proteins // F: nucleotides
# -o : Parses deflines and indexes seqIDs

formatdb -i  $INPUT_DB -p F -o T

### 1.2. Run Blast ###
# PARAMETERS:
# -p : program: "blastp", "blastn", "blastx", "tblastn", or "tblastx" 
# -d : database
# -i : query
# -e : Expectation value => -e 1e-10 [default = 0]
# -F : DUST filter => -F T/F [default = T]
# -T : Produce HTML output => -T T/F [default = F]

blastall -p tblastx -d $INPUT_DB -i $INPUT_QUERY -o $OUTPUT -T F -e $option_e -F "mS" -b 1 -v 1 -K 1 > $LOG

##############################
