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
INPUT_03="./$WORK_DIR/S03_run_blast_with_k_filter.sh"
INPUT_10="./$WORK_DIR/S04_run_blast2_with_k_filter.sh"

############################################
## 1 ## First blast run
bash $INPUT_03 WORK_DIR XX1 XX2 XX3 XX4 XX5  
############################################
## 2 ## Extract match from first blast
python S05_script_extract_match_v20_blastx.py WORK_DIR #removed $SCRIPT_PATH/05_...
############################################
## 3 ## Post processing
python S06_post_processing_of_pairwise.py WORK_DIR
python S07_format_match_get_back_nucleotides.py XX1 WORK_DIR 
############################################
## 4 ## Second blast run
bash $INPUT_10 WORK_DIR XX2 XX3 XX4 XX5
############################################
## 5 ## Extract match from second blast
python S08_script_extract_match_v20_blastx.py WORK_DIR
############################################
## 6 ## Post processing
python S09_post_processing_of_pairwise.py WORK_DIR
############################################
## 7 ## Extract RBH
python S10_compare_list_pairs_for_reciprocal_best_hits_test.py WORK_DIR
############################################
## 8 ## Post processing
python S11_post_processing_of_pairwise.py WORK_DIR
############################################
## 9 ## Conversion proteic to nucleic
python S12_prot2dna.py XX2 XX1 WORK_DIR XX3 XX4 XX5 
############################################
## 10 ## Conversion output in format zip
python S13_zip.py
