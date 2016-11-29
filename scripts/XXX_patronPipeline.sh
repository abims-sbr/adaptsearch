#$ -S /bin/bash
#$ -cwd
#$ -V

# XX1 = DB path
# XX2 = QUERY path
# XX3 = DB name
# XX4 = QUERY name
# XX5 = e_value

WORK_DIR=$1
INPUT_03="./$WORK_DIR/03_run_BLAST_with.K.filter.sh"
INPUT_10="./$WORK_DIR/10_run_BLAST2_with.K.filter.sh"

############################################
## 1 ## First blast run
bash $INPUT_03 WORK_DIR XX1 XX2 XX3 XX4 XX5  
############################################
## 2 ## Extract match from first blast
python /w/galaxy/galaxy4kevin/galaxy-dist/tools/julie/oasearch/pairwise/05_scriptExtractMatch_v20_BLASTX.py WORK_DIR
############################################
## 3 ## Post processing
python /w/galaxy/galaxy4kevin/galaxy-dist/tools/julie/oasearch/pairwise/08_1_postProcessing_of_pairwise_v3.0.py WORK_DIR
python /w/galaxy/galaxy4kevin/galaxy-dist/tools/julie/oasearch/pairwise/08_2_formatMatch_getBackNucleotides_v2.py XX1 WORK_DIR 
############################################
## 4 ## Second blast run
bash $INPUT_10 WORK_DIR XX2 XX3 XX4 XX5
############################################
## 5 ## Extract match from second blast
python /w/galaxy/galaxy4kevin/galaxy-dist/tools/julie/oasearch/pairwise/12_scriptExtractMatch_v20_BLASTX.py WORK_DIR
############################################
## 6 ## Post processing
python /w/galaxy/galaxy4kevin/galaxy-dist/tools/julie/oasearch/pairwise/14_postProcessing_of_pairwise_v3.0.py WORK_DIR
############################################
## 7 ## Extract RBH
python /w/galaxy/galaxy4kevin/galaxy-dist/tools/julie/oasearch/pairwise/16_compareListPairs_forReciprocalBestHitsTest_v3.1.py WORK_DIR
############################################
## 8 ## Post processing
python /w/galaxy/galaxy4kevin/galaxy-dist/tools/julie/oasearch/pairwise/18_postProcessing_of_pairwise_v4.0.py WORK_DIR
############################################
## 9 ## Conversion proteic to nucleic
python /w/galaxy/galaxy4kevin/galaxy-dist/tools/julie/oasearch/pairwise/24_PROT2DNA_v2.2.py XX2 XX1 WORK_DIR XX3 XX4 XX5 
############################################
## 10 ## Conversion output in format zip
python /w/galaxy/galaxy4kevin/galaxy-dist/tools/julie/oasearch/pairwise/zip.py
