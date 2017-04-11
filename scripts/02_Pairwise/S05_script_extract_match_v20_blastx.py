#!/usr/bin/env python
## AUTHOR: Eric Fontanillas
## LAST VERSION: 14/08/14 by Julie BAFFARD

### TBLASTX formatting

### MATCH = Only the first match keeped
MATCH = 0    # Only 1rst match Wanted
#MATCH = 1    # All match wanted

### SUBMATCH = several part of a same sequence match with the query
SUBMATCH = 0    # SUBMATCH NOT WANTED (only the best hit)
#SUBMATCH =1    # SUBMATCH WANTED


### NAME FORMATTING:
# [A] FORMAT QUERY NAME 1st STEP [IN DEF1]

# [B] FORMAT MATCH NAME 1st STEP [IN DEF2.1]

# [C] FORMAT MATCH NAME 2nd STEP [MIDDLE of DEF 2.3]
# [D] FORMAT QUERY NAME 2nd STEP [END of DEF 2.3]
# [E] FORMAT MATCH NAME 3rd STEP [END of DEF 2.3]

### SPECIFICITY TBLASTX (/BLASTN) formatting:
## 1/ "TBLASTX" formatting => At start of "RUN RUN RUN" change the keyword
## 2/ change line "if keyword in nextline:" in function "split_file"
## 3/ change "Strand" by "Frame" in function "get_information_on_matches"


########################################
### DEF 1. Split each "BLASTN" event ###
########################################
def split_file(path_in, keyword):

    file_in = open(path_in, "r")
    RUN = ''
    BASH1={}
    
    while 1:
        nextline = file_in.readline()

        ##################################
        ###  [A] FORMATTING QUERY NAME ###
        
        ### Get query name ###
        if nextline[0:6]=='Query=':    
            L1 = string.split(nextline, "||")
            L2 = string.split(L1[0], " ")
            query = L2[1]
            if query[-1] == "\n":
                query = query[:-1]
            
        ###  [A] END FORMATTING QUERY NAME ###
        ######################################
        
        
        ### split the file with keyword ###
        if keyword in nextline:

            # Two cases here:
            #1# If it is the first "RUN" in the block (i.e. the first occurence of "BLASTN" in the file), we have just to add the new lines in the "RUN" list ... 2nd , we have also to detect the 'key' of bash1, which is the "query" name ... and third we will have to save this "RUN" in the bash1, once we will have detected a new "RUN" (i.e. a new line beginning with "BLASTN".
            #2# If it isn't the first run, we have the save the previous "RUN" in the "bash1", before to re-initialize the RUN list (RUN =[]), before to append lines to the new "RUN"

            if RUN == '':                          # case #1#
                RUN = RUN + nextline     # we just added the first line of  the file

            else:                                # case #2# (there was a run before)
                BASH1[query] = RUN    # add the previous run to the bash
                RUN = ''                           # re-initialize the "RUN"
                RUN = RUN + nextline    # add the line starting with the keyword ("BLASTN") (except the first line of the file (the first "RUN")

        else:                          # Treatment of the subsequent lines of the one starting with the keyword ("BLASTN")  (which is not treated here but previously)
            RUN = RUN + nextline
            

        if not nextline:                                   # when no more line, we should record the last "RUN" in the bash1
            BASH1[query] = RUN                       # add the last "RUN"
            break
    
    file_in.close()
    return(BASH1)
#########################################################


################################################
### DEF2 : Parse blast output for each query ###
################################################
### detect matches (i.e. 'Sequences producing significant alignments:' ###

def detect_Matches(query, MATCH, WORK_DIR):
    F5 = open("%s/blastRun2.tmp" %WORK_DIR, 'w') 
    F5.write(bash1[query])
    F5.close()

    F6 = open("%s/blastRun2.tmp" %WORK_DIR, 'r') 
    list1 =[]
    list2 =[]

    while 1:
        nexteu = F6.readline()
        if not nexteu : break

        if "***** No hits found ******" in nexteu :
            hit = 0
            break
        
        if 'Sequences producing significant alignments:' in nexteu:
            hit = 1          
            F6.readline() # jump a line

            while 1:
                nexteu2 = F6.readline()           
                if nexteu2[0]==">": break

                ######################################
                ### [B] FORMAT MATCH NAME 1st STEP ###
               
                if nexteu2 != '\n':
                    LL1 = string.split(nexteu2, " ") # specific NORTH database names !!!!!!!
                    match = LL1[0]                     #### SOUTH databank // NORTH will have "|" separators
                    list1.append(match)

                    match2 = ">" + LL1[0]        # more complete name // still specific NORTH database names !!!!!!!
                    list2.append(match2)         #### SOUTH databank // NORTH will have "|" separators

                    if MATCH == 0:        ## Only read the 1rst line (i.e. the First Match)
                        break                
                    else:                 ## Read the other lines (i.e. All the Matches)
                        continue

                ### [B] END FORMAT MATCH NAME 1st STEP ###
                ##########################################                    
            
            break
    
    F6.close()
    return(list1, list2, hit)    # list1 = short name // list2 = more complete name
#######################################


#########################################
### DEF3 : Get Information on matches ###
#########################################
### Function used in the next function (2.3.)

def get_information_on_matches(list_of_line):   
    for line in list_of_line:

        ## Score and Expect
        if "Score" in line:
            line = line[:-1]  # remove "\n"
            S_line = string.split(line, " = ") 
            Expect = S_line[-1]                                           ## ***** Expect
            S_line2 = string.split(S_line[1], " bits ")
            Score = string.atof(S_line2[0])

        ## Identities/gaps/percent/divergence/length_matched
        elif "Identities" in line:
            line = line[:-1]  # remove "\n"
            g = 0

            if "Gaps" in line:
                pre_S_line = string.split(line, ",")
                identity_line = pre_S_line[0]
                gaps_line = pre_S_line[1]
                g = 1
            else:
                identity_line = line
                g = 0

            ## treat identity line
            S_line = string.split(identity_line, " ")

            identities = S_line[-2]                                  ## ***** identities

            S_line2 = string.split(identities, "/")
            hits = string.atof(S_line2[0])                           ## ***** hits
            length_matched = string.atof(S_line2[1])                 ## ***** length_matched
            abs_nb_differences = length_matched - hits               ## ***** abs_nb_differences


            identity_percent =   hits/length_matched * 100                   ## ***** identity_percent

            divergence_percent = abs_nb_differences/length_matched*100                  ## ***** divergence_percent

            ## treat gap line if any
            if g ==1:    # means there are gaps
                S_line3 = string.split(gaps_line, " ")
                gaps_part  = S_line3[-2]
                S_line4 = string.split(gaps_part, "/")
                gaps_number = string.atoi(S_line4[0])                                   ## ***** gaps_number
                
                real_differences = abs_nb_differences - gaps_number                         ## ***** real_differences
                real_divergence_percent = (real_differences/length_matched)*100             ## ***** real_divergence_percent
            else:
                gaps_number = 0               
                real_differences = 0
                real_divergence_percent = divergence_percent

        ## Frame
        elif "Frame" in line:
            line = line[:-1] 
            S_line = string.split(line, " = ")
            frame = S_line[1]

    list_informations=[length_matched, Expect, Score, identities, hits, identity_percent, divergence_percent,gaps_number, real_divergence_percent, frame, length_matched]

    return(list_informations)
########################################


############################
### DEF4 : get sequences ###
############################
### [+ get informations from the function 2.2.]

def get_sequences(query, list2, SUBMATCHEU,WORK_DIR):
    list_Pairwise = []

    F7 = open("%s/blastRun3.tmp" %WORK_DIR, 'w')
    F7.write(bash1[query])    # bash1[query]  ==> blast output for each query
    F7.close()
    F8 = open("%s/blastRun3.tmp" %WORK_DIR, 'r')
  
    text1 = F8.readlines()
    
    miniList = []
    for name in list2: # "list2" contains name of matched sequences (long version! the list1 is the same list but for short version names). It was previously generated by "detect_Matches" function

        l = -1
        for n in text1:
            l = l+1
            if name in n:
                i = l
                miniList.append(i)   # content positions in the list "text1", of all begining of match (e.g. >gnl|UG|Apo#S51012099 [...])

    miniList.reverse()    
    if miniList != []:
        length = len(miniList)
        ii = 0

        Listing1 = []
        while ii < length:
            iii = miniList[ii]
            entry = text1[iii:]
            text1 = text1[:iii]
            Listing1.append(entry)   # each "entry" = list of thing beginning by ">"
            ii = ii+1                # Listing1 is a table of table!!
        
        Listing1.append(text1) # "text1" = the first lines (begin with "BLASTN 2.2.1 ...]"
        Listing1.reverse()

        Listing2 = Listing1[1:]   # remove the first thing ("BLASTN ...") and keep only table beginning with a line with ">"
        SEK = len(Listing2)        
        NB_SEK = 0
        
        for e1 in Listing2:   # "Listing2" contents all the entries begining with ">"            
            NB_SEK = NB_SEK + 1
            list51 = []

            l = -1
            for line in e1:
                l = l+1
                if "Score =" in line:
                    index = l
                    list51.append(l)     # index of the lines with score

            list51.reverse()            
            Listing3 = []
            
            for i5 in list51:
                e2 = e1[i5:]
                Listing3.append(e2)
                e1 = e1[:i5]

            ######################################
            ### [C] FORMAT MATCH NAME 2nd STEP ###
            
            BigFastaName = e1      ### LIST OF LINES <=> What is remaining after removing all the hit with "Score =", so all the text comprise between ">" and the first "Score =" ==> Include Match name & "Length & empty lines
            
            SmallFastaName = BigFastaName[0] ## First line <=> MATCH NAME            
            SmallFastaName = SmallFastaName[1:-1]  ### remove ">" and "\n"

            if SmallFastaName[-1] == " ":
                SmallFastaName = SmallFastaName[:-1]
           
            PutInFastaName1 = SmallFastaName

            ### [C] END FORMAT MATCH NAME 2nd STEP ###
            ##########################################

            SUBSEK = len(Listing3)   
            NB_SUBSEK = 0
            list_inBatch = []

            ### IF NO SUBMATCH WANTED !!!! => ONLY KEEP THE FIRST HIT OF "LISTING3":
            if SUBMATCHEU == 0: # NO SUBMATCH WANTED !!!!
                Listing4 = []
                Listing4.append(Listing3[-1])   # Remove this line if submatch wanted!!!
            elif SUBMATCHEU == 1:
                Listing4 = Listing3
           
            for l in Listing4:   ## "listing3" contents               
                NB_SUBSEK = NB_SUBSEK+1
     
                ll1 = string.replace(l[0], " ", "")
                ll2 = string.replace(l[1], " ", "")
                ll3 = string.replace(l[2], " ", "")
                PutInFastaName2 = ll1[:-1] + "||" + ll2[:-1] + "||" + ll3[:-1] # match information

                seq_query = ""
                pos_query = []
                seq_match = ""
                pos_match = []

                for line in l:
                    if "Query:" in line:
                        line = string.replace(line, "    ", " ")  # remove multiple spaces in line
                        line = string.replace(line, "   ", " ")
                        line = string.replace(line, "  ", " ")
                        
                        lll1 = string.split(line, " ") # split the line, 0: "Query=", 1:start, 2:seq, 3:end

                        pos1 = lll1[1]
                        pos1 = string.atoi(pos1)
                        pos_query.append(pos1)

                        pos2 = lll1[3][:-1]
                        pos2 = string.atoi(pos2)
                        pos_query.append(pos2)
                        
                        seq = lll1[2]
                        seq_query = seq_query + seq                       

                    if "Sbjct:" in line:
                        line = string.replace(line, "    ", " ")  # remove multiple spaces in line
                        line = string.replace(line, "   ", " ")
                        line = string.replace(line, "  ", " ")
                        
                        lll2 = string.split(line, " ") # split the line, 0: "Query=", 1:start, 2:seq, 3:end

                        pos1 = lll2[1]
                        pos1 = string.atoi(pos1)
                        pos_match.append(pos1)

                        pos2 = lll2[3][:-1]
                        pos2 = string.atoi(pos2)
                        pos_match.append(pos2)
                        
                        seq = lll2[2]
                        seq_match = seq_match + seq

                ## Get the query and matched sequences and the corresponding positions
                pos_query.sort() # rank small to big
                pos_query_start = pos_query[0] # get the smaller
                pos_query_end = pos_query[-1] # get the bigger
                PutInFastaName3 = "%d...%d" %(pos_query_start, pos_query_end)

                ######################################
                ### [D] FORMAT QUERY NAME 2nd STEP ###

                FINAL_fasta_Name_Query = ">" + query + "||"+ PutInFastaName3 + "||[[%d/%d]][[%d/%d]]" %(NB_SEK, SEK, NB_SUBSEK,SUBSEK)

                ### [D] END FORMAT QUERY NAME 2nd STEP ###
                ##########################################

                pos_match.sort()
                pos_match_start = pos_match[0]
                pos_match_end = pos_match[-1]
                PutInFastaName4 = "%d...%d" %(pos_match_start, pos_match_end)

                ######################################
                ### [E] FORMAT MATCH NAME 3rd STEP ###
               
                FINAL_fasta_Name_Match = ">" + PutInFastaName1 + "||" + PutInFastaName4 + "||[[%d/%d]][[%d/%d]]" %(NB_SEK, SEK, NB_SUBSEK,SUBSEK)

                ### [E] END FORMAT MATCH NAME 3rd STEP ###
                ##########################################

                Pairwise = [FINAL_fasta_Name_Query , seq_query , FINAL_fasta_Name_Match , seq_match]  # list with 4 members
                list_Pairwise.append(Pairwise)

                ### Get informations about matches
                list_info = get_information_on_matches(l)    ### DEF3 ###
                                    
    F8.close()
    return(list_Pairwise, list_info)
#########################################


######################
### 2. RUN RUN RUN ###
######################
import string, os, time, re, sys

## 1 ## INPUT/OUTPUT
SHORT_FILE = sys.argv[1] #short-name-query_short-name-db

path_in = "%s/04_outputBlast_%s.txt" %(SHORT_FILE, SHORT_FILE) 
file_out = open("%s/06_PairwiseMatch_%s.fasta" %(SHORT_FILE, SHORT_FILE),"w")

## 2 ## RUN
## create Bash1 ##
bash1 = split_file(path_in, "TBLASTX") ## DEF1 ##

## detect and save match ##
list_hits =[]
list_no_hits = []
j= 0
k = 0
lene = len(bash1.keys())
for query in bash1.keys():
    j  = j+1
    
    ## 2.1. detect matches ##
    list_match, list_match2, hit=detect_Matches(query, MATCH, SHORT_FILE)    ### DEF2 ###
    
    if hit == 1:   # match(es)
        list_hits.append(query)
    if hit == 0: # no match for that sequence
        list_no_hits.append(query)
    
    ## 2.2. get sequences ##
    if hit ==1:
        list_pairwiseMatch, list_info = get_sequences(query, list_match2, SUBMATCH, SHORT_FILE)       ### FUNCTION ###

        # divergencve
        divergence = list_info[6]
        # gap number
        gap_number = list_info[7]
        # real divergence (divergence without accounting INDELs)
        real_divergence = list_info[8]
        # length matched
        length_matched = list_info[10]

        ### WRITE PAIRWISE ALIGNMENT IN OUTPUT FILES
        for pairwise in list_pairwiseMatch:
            k = k+1

            query_name = pairwise[0]
            query_seq = pairwise[1]
            match_name = pairwise[2]
            match_seq = pairwise[3]

            len_query_seq = len(query_seq)

            Lis1 = string.split(query_name, "||")
            short_query_name = Lis1[0]
            Lis2 = string.split(match_name, "||")
            short_match_name = Lis2[0]

            # If NO CONTROL FOR LENGTH, USE THE FOLLOWING LINES INSTEAD:
            
            file_out.write("%s||%s||%s||%s||%s" %(query_name,divergence,gap_number,real_divergence,length_matched))
            file_out.write("\n")
            file_out.write("%s" %query_seq)
            file_out.write("\n")

            file_out.write("%s||%s||%s||%s||%s" %(match_name,divergence,gap_number,real_divergence,length_matched))
            file_out.write("\n")
            file_out.write("%s" %match_seq)
            file_out.write("\n")
    
file_out.close()
