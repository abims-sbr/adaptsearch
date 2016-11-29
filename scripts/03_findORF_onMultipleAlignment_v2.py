#!/usr/bin/python
## Author: Eric Fontanillas
## Last modification: 13/10/2010
## Subject: Predict potential ORF on the basis of 2 criteria + 1 optional criteria
            ## CRITERIA 1 ## Longest part of the alignment of sequence without codon stop "*", tested in the 3 potential ORF
            ## CRITERIA 2 ## This longest part should be > 150nc or 50aa
            ## CRITERIA 3 ## [OPTIONNAL] A codon start "M" should be present in this longuest part, before the last 50 aa
                                 ## OUTPUTs "05_CDS_aa" & "05_CDS_nuc" => NOT INCLUDE THIS CRITERIA
                                 ## OUTPUTs "06_CDS_with_M_aa" & "06_CDS_with_M_nuc" => INCLUDE THIS CRITERIA


###############################
##### DEF 0 : Dico fasta  #####
###############################
def dico(fasta_file_path):
    F2 = open(fasta_file_path, "r")
    dicoco = {}
    while 1:
        next2 = F2.readline()
        if not next2:
            break
        if next2[0] == ">":
            fasta_name_query = next2[:-1]
            Sn = string.split(fasta_name_query, "||")
            fasta_name_query = Sn[0]
            next3 = F2.readline()
            fasta_seq_query = next3[:-1]
            dicoco[fasta_name_query]=fasta_seq_query
    F2.close()
    return(dicoco)
###################################################################################


#####################
###### DEF 1 ########
#####################
#################################
###  Create bash for genetic code
###       KEY = codon
###       VALUE = Amino Acid
#################################
def code_universel(F1):
    bash_codeUniversel = {}
    while 1:
        next = F1.readline()
        if not next: break
        L1 = string.split(next, " ")
        length1 = len(L1)
        if length1 == 3:
            key = L1[0]
            value = L1[2][:-1]
            bash_codeUniversel[key] = value
        else:
            key =  L1[0]
            value = L1[2]
            bash_codeUniversel[key] = value
    #print bash_codeUniversel
    F1.close()
    return(bash_codeUniversel)
###########################################################



#################
##### DEF 2 #####
#################
############################################################
### Test if the sequence is a multiple of 3, and if not correct the sequence to become a multiple of 3 ###
### !!!!!!!!!!!!!!!!!!!!! WEAKNESS OF THAT APPROACH = I remove extra base(s) at the end of the sequence ==> I can lost a codon, when I test ORF (as I will decay the ORF)
############################################################
def multiple3(seq):
    leng = len(seq)
    #print "\nINITIAL LENGTH = %d" %leng
    modulo = leng%3
    if modulo == 0:   # the results of dividing leng per 3 is an integer
        new_seq = seq        
    elif modulo == 1:    # means 1 extra nc (nucleotid) needs to be removed (the remaining of modulo indicate the part which is non-dividable per 3)
        new_seq = seq[:-1]  # remove the last nc
    elif modulo == 2:  # means 2 extra nc (nucleotid) needs to be removed (the remaining of modulo indicate the part which is non-dividable per 3)
        new_seq = seq[:-2]  # remove the 2 last nc
    len1 = len(new_seq)
    #print "NEW LENGTH = %d\n" %len1    
    return(new_seq, modulo)
##########################################################



#####################
####### DEF 3 #######
#####################
## GET ORF:
##         - MONO SEQUENCE BASED : Based on 1 sequence only!
##         - CRITERIA: Minimum number of codon stop in the ORF


##### DEF 3 - PART 1 -  ##### (Require ### DEF 3 - PART 2 - ###)   ==> WARNING ===> UN-TESTED FUNCTION !!!!!!!!!!!!!!!!!
###############################
def get_ORF_codonStopNumberCriteria(seq, bash_codeUniversel):    
    seq_dna = ""
    seq_aa = ""
    codon_stop_nb=0
    i = 0
    len1 = len(seq)
    while i < len1:
        base1 = seq[i]
        base1 = string.capitalize(base1)
        base2 = seq[i+1]
        base2 = string.capitalize(base2)
        base3 = seq[i+2]
        base3 = string.capitalize(base3)
        
        codon = base1+base2+base3
        seq_dna = seq_dna + codon
        codon = string.replace(codon, "T", "U")

        if codon in bash_codeUniversel.keys():
            aa = bash_codeUniversel[codon]
            seq_aa = seq_aa + aa            
            if aa == "*":
                codon_stop_nb = codon_stop_nb+1
                #print codon
        else:
            seq_aa = seq_aa +"?"    ### Take account for gap "-" and "N"
        i = i + 3
    
    return(seq_dna, seq_aa, codon_stop_nb)

##### DEF 3 - PART 2 - #####
##############################
## Tests the 3 ORF and returns a list of ORF OK. WARNING = SEVERAL ORF MAYBE OK!!!!
def find_good_ORF_criteria_1(new_seq, bash_codeUniversel):  
    len1 = len(new_seq)
    
    ####################
    ### 1. TEST ORF1 ###
    ####################
    print "TEST ORF1"
    seq = new_seq  # no change
    seq_dna_ORF1, seq_aa_ORF1, codon_stop_nb_ORF1 = get_ORF(seq, bash_codeUniversel)       ### DEF 3 - PART 1 - ###
    
    ####################
    ### 2. TEST ORF2 ###
    ####################
    len2 = len1
    print "TEST ORF2"
    seq = new_seq[1:-2]  # remove 1 position (nc) at start and 2 positions at end (to respect codon structure)
    len2 = len2 - 3 # <=> removing 1 codon
    seq_dna_ORF2, seq_aa_ORF2, codon_stop_nb_ORF2 = get_ORF(seq, bash_codeUniversel)       ### DEF 3 - PART 1 - ###

    ####################
    ### 3. TEST ORF3 ###
    ####################
    len3 = len1
    print "TEST ORF3"
    seq = new_seq[2:-1]  # remove 2 positions (nc) at start and 1 position at end (to respect codon structure)
    len3 = len3 - 3   # <=> removing 1 codon

    seq_dna_ORF3, seq_aa_ORF3, codon_stop_nb_ORF3 = get_ORF(seq, bash_codeUniversel)       ### DEF 3 - PART 1 - ###
    
    #####################################################################################
    ### 4. list Best ORF (on the basis of the nb of codon stop detected for each ORF) ###
    #####################################################################################
    ## !!!! SEVERAL ORF MAYBE GOOD (AND NOT ONLY 1) !!!
    list_ORFok = [codon_stop_nb_ORF1,codon_stop_nb_ORF2,codon_stop_nb_ORF3]

    
    ################
    ### 6.RETURN ###
    ################       
    return(list_ORFok, seq_dna_ORF1 , seq_aa_ORF1, seq_dna_ORF2 , seq_aa_ORF2, seq_dna_ORF3 , seq_aa_ORF3)
##########################################################




###################
###### DEF 4 ######
###################
## GET ORF: (Require ### DEF 3.2. - PART 2 - ###)
##         - MONO SEQUENCE BASED : Based on 1 sequence only!
##         - CRITERIA1: Longuest part of the sequence without Codon Stop
##         - CRITERIA2: That longuest part should be at least longer than 50 aa (150bp) ==> if not: probably not a coding sequence, or a too truncated coding sequence

##         - CRITERIA3: [OPTIONNAL] A codon start "M" should be present in this longuest part, before the last 50 aa
                                 ## OUTPUTs "05_CDS_aa" & "05_CDS_nuc" => NOT INCLUDE THIS CRITERIA
                                 ## OUTPUTs "06_CDS_with_M_aa" & "06_CDS_with_M_nuc" => INCLUDE THIS CRITERIA


###### DEF 4 - Part 1 - ###### (Require ### DEF 4. - PART 2 - ###)   ==> WARNING ===> UN-TESTED FUNCTION !!!!!!!!!!!!!!!!!
##############################
def get_ORF_LonguestSeqWithoutCodonStopCriteria(seq, bash_codeUniversel):
    seq_dna = ""
    seq_aa = ""
    i = 0
    len1 = len(seq)
    while i < len1:
        base1 = seq[i]
        base1 = string.capitalize(base1)
        base2 = seq[i+1]
        base2 = string.capitalize(base2)
        base3 = seq[i+2]
        base3 = string.capitalize(base3)
        
        codon = base1+base2+base3
        seq_dna = seq_dna + codon
        #codon = string.replace(codon, "-", "N")
        codon = string.replace(codon, "T", "U")

        if codon in bash_codeUniversel.keys():
            aa = bash_codeUniversel[codon]
            seq_aa = seq_aa + aa            
        else:
            seq_aa = seq_aa +"?"    ### Take account for gap "-" and "N"
        i = i + 3

    ## test the length of the longuest part of the sequence without "*" (i.e. codon stop)
    if "*" in seq_aa:
        S1 = string.split(seq_aa, "*")
        MAX_AA_LENGTH=0
        LONGUEST_PART_WHITHOUT_CODON_STOP = ""
        for subsequence in S1:
            if len(subsequence) > MAX_AA_LENGTH:
                LONGUEST_PART_WHITHOUT_CODON_STOP = subsequence
                MAX_AA_LENGTH = len(subsequence)
    else:
        MAX_AA_LENGTH = len(seq_aa)
        LONGUEST_PART_WHITHOUT_CODON_STOP = seq_aa
        
    
    return(seq_dna, seq_aa, LONGUEST_PART_WHITHOUT_CODON_STOP, MAX_AA_LENGTH)
    

##### DEF 4 - Part 2 - #####
############################
## Tests the 3 ORF and returns a list of ORF OK. WARNING = SEVERAL ORF MAYBE OK!!!!
## Codon structure is preserved when removing 1 nc at start (so 2 nc removed at the end) or when removing 2 nc at start (so 1 nc removed at the end)
def find_good_ORF_criteria_2(new_seq, bash_codeUniversel):
    len1 = len(new_seq)
    ####################
    ### 1. TEST ORF1 ###
    ####################
    print "TEST ORF1"
    seq = new_seq  # no change
    seq_dna_ORF1, seq_aa_ORF1, LONGUEST_PART_WHITHOUT_CODON_STOP_ORF1, MAX_AA_LENGTH_ORF1 = get_ORF(seq, bash_codeUniversel)       ### DEF 4 - Part 1 - ###
    
    ####################
    ### 2. TEST ORF2 ###
    ####################
    len2 = len1
    print "TEST ORF2"
    seq = new_seq[1:-2]  # remove 1 position (nc) at start and 2 positions at end (to respect codon structure)
    len2 = len2 - 3 # <=> removing 1 codon
    seq_dna_ORF2, seq_aa_ORF2, LONGUEST_PART_WHITHOUT_CODON_STOP_ORF2, MAX_AA_LENGTH_ORF2 = get_ORF(seq, bash_codeUniversel)       ### DEF 4 - Part 1 - ###

    ####################
    ### 3. TEST ORF3 ###
    ####################
    len3 = len1
    print "TEST ORF3"
    seq = new_seq[2:-1]  # remove 2 positions (nc) at start and 1 position at end (to respect codon structure)
    len3 = len3 - 3   # <=> removing 1 codon
    seq_dna_ORF3, seq_aa_ORF3, LONGUEST_PART_WHITHOUT_CODON_STOP_ORF3, MAX_AA_LENGTH_ORF3 = get_ORF(seq, bash_codeUniversel)       ### DEF 4 - Part 1 - ###
    
    #####################################################################################
    ### 4. list Best ORF (on the basis of the nb of codon stop detected for each ORF) ###
    #####################################################################################
    ## !!!! SEVERAL ORF MAYBE GOOD (AND NOT ONLY 1) !!!
    list_ORFok = [MAX_AA_LENGTH_ORF1,MAX_AA_LENGTH_ORF2,MAX_AA_LENGTH_ORF3]
    
    ################
    ### 6.RETURN ###
    ################       
    return(list_ORFok, seq_dna_ORF1, seq_aa_ORF1, seq_dna_ORF2, seq_aa_ORF2, seq_dna_ORF3, seq_aa_ORF3)
##########################################################


###################
###### DEF 5 ######
###################
## GET ORF: (Require ### DEF 3.2. - PART 2 - ###)
##         - MULTIPLE SEQUENCE BASED : Based on ALIGNMENT of several sequences
##         - CRITERIA1: Get the segment in the alignment with no codon stop



###### DEF 5 - Part 1 - ######
##############################
def simply_get_ORF(seq_dna, bash_codeUniversel):
    seq_aa = ""
    i = 0
    len1 = len(seq_dna)
    while i < len1:
        base1 = seq_dna[i]
        base1 = string.capitalize(base1)
        base2 = seq_dna[i+1]
        base2 = string.capitalize(base2)
        base3 = seq_dna[i+2]
        base3 = string.capitalize(base3)
        
        codon = base1+base2+base3
        codon = string.replace(codon, "T", "U")

        if codon in bash_codeUniversel.keys():
            aa = bash_codeUniversel[codon]
            seq_aa = seq_aa + aa            
        else:
            seq_aa = seq_aa +"?"    ### Take account for gap "-" and "N"
        i = i + 3

    return(seq_aa)


###### DEF 5 - Part 2 - ######
##############################
def find_good_ORF_criteria_3(bash_aligned_nc_seq, bash_codeUniversel):  
    
    ## 1 ## Get the list of aligned aa seq for the 3 ORF:
    bash_of_aligned_aa_seq_3ORF = {}
    bash_of_aligned_nuc_seq_3ORF = {}
    BEST_LONGUEST_SUBSEQUENCE_LIST_POSITION = []
    for fasta_name in bash_aligned_nc_seq.keys():
        ## 1.1. ## Get the raw sequence
        sequence_nc = bash_aligned_nc_seq[fasta_name]

        ## 1.2. ## Check whether the sequence is multiple of 3, and correct it if not:
        new_sequence_nc, modulo = multiple3(sequence_nc)  ### DEF 2 ###
        
        ## 1.3. ## Get the 3 ORFs (nuc) for each sequence
        seq_nuc_ORF1 = new_sequence_nc
        seq_nuc_ORF2 = new_sequence_nc[1:-2]
        seq_nuc_ORF3 = new_sequence_nc[2:-1]
        
        LIST_3_ORF_nuc = [seq_nuc_ORF1, seq_nuc_ORF2, seq_nuc_ORF3]
        bash_of_aligned_nuc_seq_3ORF[fasta_name] = LIST_3_ORF_nuc     ### For each seq of the multialignment => give the 3 ORF (in nuc)
        

        ## 1.4. ## Get the 3 ORFs (aa) for each sequence
        seq_prot_ORF1 = simply_get_ORF(seq_nuc_ORF1,bash_codeUniversel)   ### DEF 5 - Part 1 - ##
        seq_prot_ORF2 = simply_get_ORF(seq_nuc_ORF2,bash_codeUniversel)   ### DEF 5 - Part 1 - ##
        seq_prot_ORF3 = simply_get_ORF(seq_nuc_ORF3,bash_codeUniversel)   ### DEF 5 - Part 1 - ##

        LIST_3_ORF_aa = [seq_prot_ORF1, seq_prot_ORF2, seq_prot_ORF3]
        bash_of_aligned_aa_seq_3ORF[fasta_name] = LIST_3_ORF_aa     ### For each seq of the multialignment => give the 3 ORFs (in aa)
        
        

    ## 2 ## Test for the best ORF (Get the longuest segment in the alignment with no codon stop ... for each ORF ... the longuest should give the ORF)
    BEST_MAX = 0
    for i in [0,1,2]:   ### Test the 3 ORFs 
        ORF_Aligned_aa = []
        ORF_Aligned_nuc = []
                

        ## 2.1 ## Get the alignment of sequence for a given ORF
        ## Compare the 1rst ORF between all sequence => list them in ORF_Aligned_aa // them do the same for the second ORF, and them the 3rd
        for fasta_name in bash_of_aligned_aa_seq_3ORF.keys(): 
            ORFsequence = bash_of_aligned_aa_seq_3ORF[fasta_name][i]
            aa_length = len(ORFsequence)
            ORF_Aligned_aa.append(ORFsequence)   ### List of all sequences in the ORF nb "i" =

        n = i+1
#         print "ORF %d = " %n
#         print ORF_Aligned_aa
#         for sek in ORF_Aligned_aa:
#             print sek
#         print len(ORF_Aligned_aa)

        
        for fasta_name in bash_of_aligned_nuc_seq_3ORF.keys(): 
            ORFsequence = bash_of_aligned_nuc_seq_3ORF[fasta_name][i]
            nuc_length = len(ORFsequence)
            ORF_Aligned_nuc.append(ORFsequence)   ### List of all sequences in the ORF nb "i" =

#         print "ORF %d = " %n
#         print ORF_Aligned_nuc
#         for sek in ORF_Aligned_nuc:
#             print sek
#         print len(ORF_Aligned_nuc)


        ## 2.2 ## Get the list of sublist of positions whithout codon stop in the alignment
        ## For each ORF, now we have the list of sequences available (i.e. THE ALIGNMENT IN A GIVEN ORF)
        ## Next step is to get the longuest subsequence whithout stop
        ## We will explore the presence of stop "*" in each column of the alignment, and get the positions of the segments between the positions with "*"
        MAX_LENGTH = 0
        LONGUEST_SEGMENT_UNSTOPPED = ""        
        j = 0 # Start from first position in alignment       
        List_of_List_subsequences = []
        List_positions_subsequence = []
        while j < aa_length:                
                column = []
                for seq in ORF_Aligned_aa:
                    column.append(seq[j])
                j = j+1
                if "*" in column:
                    #print "STOP FOUND"
                    List_of_List_subsequences.append(List_positions_subsequence) ## Add previous list of positions
                    List_positions_subsequence = []                              ## Re-initialyse list of positions
                else:
                    List_positions_subsequence.append(j)
        #print List_of_List_subsequences

        
        ## 2.3 ## Among all the sublists (separated by column with codon stop "*"), get the longuest one (BETTER SEGMENT for a given ORF)
        LONGUEST_SUBSEQUENCE_LIST_POSITION = []
        MAX=0
        for sublist in List_of_List_subsequences:
            if len(sublist) > MAX and len(sublist) > MINIMAL_CDS_LENGTH:
                MAX = len(sublist)
                LONGUEST_SUBSEQUENCE_LIST_POSITION = sublist
        
        ## 2.4. ## Test if the longuest subsequence start exactly at the beginning of the original sequence (i.e. means the ORF maybe truncated)
        if LONGUEST_SUBSEQUENCE_LIST_POSITION != []:
            print LONGUEST_SUBSEQUENCE_LIST_POSITION[0]
            if LONGUEST_SUBSEQUENCE_LIST_POSITION[0] == 0:
                CDS_maybe_truncated = 1
            else:
                CDS_maybe_truncated = 0
        else:
            CDS_maybe_truncated = 0
            

        ## 2.5 ## Test if this BETTER SEGMENT for a given ORF, is the better than the one for the other ORF (GET THE BEST ORF)
        ## Test whether it is the better ORF
        if MAX > BEST_MAX:
            BEST_MAX = MAX
            BEST_ORF = i+1
            BEST_LONGUEST_SUBSEQUENCE_LIST_POSITION = LONGUEST_SUBSEQUENCE_LIST_POSITION

    #print "The best ORF is the nb %d" %BEST_ORF
    #print "The longuest segment whithout stop in this ORF is the nb %s" %BEST_MAX


    ## 3 ## ONCE we have this better segment (BEST CODING SEGMENT)
    ## ==> GET THE STARTING and ENDING POSITIONS (in aa position and in nuc position)
    ## And get the INDEX of the best ORF [0, 1, or 2]
    if BEST_LONGUEST_SUBSEQUENCE_LIST_POSITION != []:
        pos_MIN_aa = BEST_LONGUEST_SUBSEQUENCE_LIST_POSITION[0]
        pos_MIN_aa = pos_MIN_aa - 1
        pos_MAX_aa = BEST_LONGUEST_SUBSEQUENCE_LIST_POSITION[-1]


        BESTORF_bash_of_aligned_aa_seq = {}
        BESTORF_bash_of_aligned_aa_seq_CODING = {}
        for fasta_name in bash_of_aligned_aa_seq_3ORF.keys():
            index_BEST_ORF = BEST_ORF-1  ### cause list going from 0 to 2 in LIST_3_ORF, while the ORF nb is indexed from 1 to 3
            seq = bash_of_aligned_aa_seq_3ORF[fasta_name][index_BEST_ORF]
            seq_coding = seq[pos_MIN_aa:pos_MAX_aa]
            #print seq_coding
            BESTORF_bash_of_aligned_aa_seq[fasta_name] = seq
            BESTORF_bash_of_aligned_aa_seq_CODING[fasta_name] = seq_coding

        ## 4 ## Get the corresponding position (START/END of BEST CODING SEGMENT) for nucleotides alignment
        pos_MIN_nuc = pos_MIN_aa * 3
        pos_MAX_nuc = pos_MAX_aa * 3

        BESTORF_bash_aligned_nc_seq = {}
        BESTORF_bash_aligned_nc_seq_CODING = {}    
        for fasta_name in bash_aligned_nc_seq.keys():
            seq = bash_of_aligned_nuc_seq_3ORF[fasta_name][index_BEST_ORF]
            seq_coding = seq[pos_MIN_nuc:pos_MAX_nuc]
            #print seq_coding
            BESTORF_bash_aligned_nc_seq[fasta_name] = seq
            BESTORF_bash_aligned_nc_seq_CODING[fasta_name] = seq_coding

    else:
        print "NO BEST CDS FOUND!!!!"
        BESTORF_bash_aligned_nc_seq = {}
        BESTORF_bash_aligned_nc_seq_CODING = {}
        BESTORF_bash_of_aligned_aa_seq = {}
        BESTORF_bash_of_aligned_aa_seq_CODING ={}



    ### Check whether their is a "M" or not, and if at least 1 "M" is present, that it is not in the last 50 aa
    ###########################################################################################################

    BESTORF_bash_of_aligned_aa_seq_CDS_with_M = {}
    BESTORF_bash_of_aligned_nuc_seq_CDS_with_M = {}
    
    Ortho = 0
    for fasta_name in BESTORF_bash_of_aligned_aa_seq_CODING.keys():
        seq_aa = BESTORF_bash_of_aligned_aa_seq_CODING[fasta_name]
        Ortho = detect_Methionine(seq_aa, Ortho)   ### DEF7 ###

    ## CASE 1: A "M" is present and correctly localized (not in last 50 aa)
    if Ortho == 1:
        BESTORF_bash_of_aligned_aa_seq_CDS_with_M = BESTORF_bash_of_aligned_aa_seq_CODING
        BESTORF_bash_of_aligned_nuc_seq_CDS_with_M = BESTORF_bash_aligned_nc_seq_CODING

    ## CASE 2: in case the CDS is truncated, so the "M" is maybe missing:
    if Ortho == 0 and CDS_maybe_truncated == 1: 
        BESTORF_bash_of_aligned_aa_seq_CDS_with_M = BESTORF_bash_of_aligned_aa_seq_CODING
        BESTORF_bash_of_aligned_nuc_seq_CDS_with_M = BESTORF_bash_aligned_nc_seq_CODING

    ## CASE 3: CDS not truncated AND no "M" found in good position (i.e. before the last 50 aa):
        ## => the 2 bash "CDS_with_M" are left empty ("{}")
    
    return(BESTORF_bash_aligned_nc_seq,  BESTORF_bash_aligned_nc_seq_CODING, BESTORF_bash_of_aligned_nuc_seq_CDS_with_M, BESTORF_bash_of_aligned_aa_seq, BESTORF_bash_of_aligned_aa_seq_CODING, BESTORF_bash_of_aligned_aa_seq_CDS_with_M)    
##########################################################



###################
###### DEF 6 ######
###################
## Detect all indices corresponding to all occurance of a substring in a string
def allindices(string, sub):
    listindex=[]
    offset=0
    i = string.find(sub, offset)
    while i >= 0:
        listindex.append(i)
        i = string.find(sub, i + 1)
    return listindex
######################################################

###################
###### DEF 7 ######
###################
## Detect if methionin in the aa sequence
def detect_Methionine(seq_aa, Ortho):

    ln = len(seq_aa)

    CUTOFF_Last_50aa = ln -50

    #Ortho = 0  ## means orthologs not found
    
    ## Find all indices of occurances of "M" in a string of aa
    list_indices = allindices(seq_aa, "M")  ### DEF6 ###
         
    ## If some "M" are present, find whether the first "M" found is not in the 50 last aa (indice < CUTOFF_Last_50aa) ==> in this case: maybenot a CDS
    if list_indices != []:
        first_M = list_indices[0]
        print first_M
        print CUTOFF_Last_50aa
        
        if first_M < CUTOFF_Last_50aa:
            Ortho = 1  ## means orthologs found

    return(Ortho)
####################################


#######################
##### RUN RUN RUN #####
#######################
import string, os, time, re

### 0 ### PARAMETERS
MINIMAL_CDS_LENGTH = 50  ## in aa number

### 1 ### OPEN FILES

## INPUT
Path_IN = "../05_Align_Orthologs_with_BLASTALIGN/28_Alignments/"
#Path_IN = "01_test/"
L_IN = os.listdir(Path_IN)

F2 = open("02_input_code_universel_modified.txt", 'r')

F3 = open("07_files_with_no_CDS.txt", "w")

LOG = open("07_FilesTreated_CDSfound.log", "w")

## OUTPUT (open and clean it!)
Path_OUT1 = "04_BEST_ORF_nuc"
Path_OUT2 = "04_BEST_ORF_aa"

Path_OUT3 = "05_CDS_nuc"
Path_OUT4 = "05_CDS_aa"

Path_OUT5 = "06_CDS_with_M_nuc"
Path_OUT6 = "06_CDS_with_M_aa"

os.system("rm 04_BEST_ORF_nuc/*")
os.system("rm 04_BEST_ORF_aa/*")
os.system("rm 05_CDS_nuc/*")
os.system("rm 05_CDS_aa/*")
os.system("rm 06_CDS_with_M_nuc/*")
os.system("rm 06_CDS_with_M_aa/*")


### 2 ### Get Universal Code
bash_codeUniversel = code_universel(F2)  ### DEF1 ###
print bash_codeUniversel
F2.close()

### 3 ### Get the Bash corresponding to an alignment file in fasta format

count_file_processed = 0
count_file_with_CDS = 0
count_file_without_CDS = 0
count_file_with_CDS_plus_M = 0

for file in L_IN:
    count_file_processed = count_file_processed + 1
    print file
    
    fasta_file_path = "%s/%s" %(Path_IN, file)

    bash_fasta = dico(fasta_file_path)   ### DEF 0 ###

    #print bash_fasta.keys()

    BESTORF_nuc, BESTORF_nuc_CODING, BESTORF_nuc_CDS_with_M, BESTORF_aa, BESTORF_aa_CODING, BESTORF_aa_CDS_with_M  = find_good_ORF_criteria_3(bash_fasta, bash_codeUniversel)   ### DEF 5 - PART 2 - ###
    
    print "ORF detection done ..."


    ## a ## OUTPUT BESTORF_nuc
    if BESTORF_nuc != {}:
        count_file_with_CDS = count_file_with_CDS +1
        OUT1 = open("%s/%s" %(Path_OUT1,file), "w")
        for fasta_name in BESTORF_nuc.keys():
            seq = BESTORF_nuc[fasta_name]
            OUT1.write("%s\n" %fasta_name)
            OUT1.write("%s\n" %seq)
        OUT1.close()
    else:
        count_file_without_CDS = count_file_without_CDS + 1
        F3.write("%s\n" %file)
        
    
    ## b ## OUTPUT BESTORF_nuc_CODING  ===> THE MOST INTERESTING!!!
    if BESTORF_aa != {}:
        #count_file_with_CDS = count_file_with_CDS +1
        OUT2 = open("%s/%s" %(Path_OUT2,file), "w")
        for fasta_name in BESTORF_aa.keys():
            seq = BESTORF_aa[fasta_name]
            OUT2.write("%s\n" %fasta_name)
            OUT2.write("%s\n" %seq)
        OUT2.close()    

    ## c ## OUTPUT BESTORF_aa
    if BESTORF_nuc_CODING != {}:
        #count_file_with_CDS = count_file_with_CDS +1
        OUT3 = open("%s/%s" %(Path_OUT3,file), "w")
        for fasta_name in BESTORF_nuc_CODING.keys():
            seq = BESTORF_nuc_CODING[fasta_name]
            OUT3.write("%s\n" %fasta_name)
            OUT3.write("%s\n" %seq)
        OUT3.close()

    ## d ## OUTPUT BESTORF_aa_CODING
    if BESTORF_aa_CODING != {}:
        #count_file_with_CDS = count_file_with_CDS +1
        OUT4 = open("%s/%s" %(Path_OUT4,file), "w")
        for fasta_name in BESTORF_aa_CODING.keys():
            seq = BESTORF_aa_CODING[fasta_name]
            OUT4.write("%s\n" %fasta_name)
            OUT4.write("%s\n" %seq)
        OUT4.close()    

    ## e ## OUTPUT BESTORF_nuc_CDS_with_M
    if BESTORF_nuc_CDS_with_M != {}:
        count_file_with_CDS_plus_M = count_file_with_CDS_plus_M + 1
        OUT5 = open("%s/%s" %(Path_OUT5,file), "w")
        for fasta_name in BESTORF_nuc_CDS_with_M.keys():
            seq = BESTORF_nuc_CDS_with_M[fasta_name]
            OUT5.write("%s\n" %fasta_name)
            OUT5.write("%s\n" %seq)
        OUT5.close()    

    ## f ## OUTPUT BESTORF_aa_CDS_with_M
    if BESTORF_aa_CDS_with_M != {}:
        OUT6 = open("%s/%s" %(Path_OUT6,file), "w")
        for fasta_name in BESTORF_aa_CDS_with_M.keys():
            seq = BESTORF_aa_CDS_with_M[fasta_name]
            OUT6.write("%s\n" %fasta_name)
            OUT6.write("%s\n" %seq)
        OUT6.close()    


LOG.write("\n\nSUMARIZE CDS DETECTION:\n")
LOG.write("\tFiles processed: %d\n" %count_file_processed)
LOG.write("\tFiles with CDS: %d\n" %count_file_with_CDS)
LOG.write("\tFiles with CDS plus M (codon start): %d\n" %count_file_with_CDS_plus_M)
LOG.write("\tFiles without CDS: %d\n" %count_file_without_CDS)


F3.close()
