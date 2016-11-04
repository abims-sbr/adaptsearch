#!/bin/bash

# XX1 = DB path
# XX2 = QUERY path
# XX3 = WORK_DIR


############################################
## 1 ## First blast run
./03_run_BLAST_with.K.filter.sh XX1 XX2
############################################
## 2 ## Extract match from first blast 
./05_scriptExtractMatch_v20_BLASTX.py "."
############################################
## 3 ## Post processing
./08_1_postProcessing_of_pairwise_v3.0.py "."
./08_2_formatMatch_getBackNucleotides_v2.py XX1 "."
############################################

./10_run_BLAST2_with.K.filter.sh XX2

./12_scriptExtractMatch_v20_BLASTX.py "."

./14_postProcessing_of_pairwise_v3.0.py "."

./16_compareListPairs_forReciprocalBestHitsTest_v3.1.py "."

./18_postProcessing_of_pairwise_v4.0.py "."

./20_get_divergence_per_pairwise_v2.py "."

./24_PROT2DNA_v2.2.py XX2 XX1 "."
