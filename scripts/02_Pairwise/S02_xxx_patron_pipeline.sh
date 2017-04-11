#$ -S /bin/env bash
#$ -cwd
#$ -V

# XX1 = DB path
# XX2 = QUERY path
# XX3 = DB name
# XX4 = QUERY name
# XX5 = e_value

WORK_DIR=$1
#SCRIPT_PATH=/home/umr7144/abice/vmataigne/Documents/AdaptSearch/adaptsearch-master/scripts
INPUT_03="./$WORK_DIR/T2S03_run_blast_with_k_filter.sh"
INPUT_10="./$WORK_DIR/T2S04_run_blast2_with_k_filter.sh"

############################################
## 1 ## First blast run
bash $INPUT_03 WORK_DIR XX1 XX2 XX3 XX4 XX5  
############################################
## 2 ## Extract match from first blast
python T2S05_script_extract_match_v20_blastx.py WORK_DIR #removed $SCRIPT_PATH/05_...
############################################
## 3 ## Post processing
python T2S06_post_processing_of_pairwise_v3.0.py WORK_DIR
python T2S07_format_match_get_back_nucleotides_v2.py XX1 WORK_DIR 
############################################
## 4 ## Second blast run
bash $INPUT_10 WORK_DIR XX2 XX3 XX4 XX5
############################################
## 5 ## Extract match from second blast
python T2S08_script_extract_match_v20_blastx.py WORK_DIR
############################################
## 6 ## Post processing
python T2S09_post_processing_of_pairwise_v3.0.py WORK_DIR
############################################
## 7 ## Extract RBH
python T2S10_compare_list_pairs_for_reciprocal_best_hits_test_v3.1.py WORK_DIR
############################################
## 8 ## Post processing
python T2S11_post_processing_of_pairwise_v4.0.py WORK_DIR
############################################
## 9 ## Conversion proteic to nucleic
python T2S12_prot2dna_v2.2.py XX2 XX1 WORK_DIR XX3 XX4 XX5 
############################################
## 10 ## Conversion output in format zip
python T2S13_zip.py
