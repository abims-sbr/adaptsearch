import string

# Used in S05 and
def split_file(path_in, keyword):

    file_in = open(path_in, "r")    
    RUN = ''
    BASH1={}
    
    with open(path_in, "r") as file_in:    
        for nextline in file_in.readlines():       

            ##################################
            ###  [A] FORMATTING QUERY NAME ###

            # Get query name
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
            

        if RUN:                           
            BASH1[query] = RUN                       # add the last "RUN"
   
    
    return(BASH1)

def detect_Matches(query, MATCH, WORK_DIR, bash1):
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

def get_information_on_matches(list_of_line):

    for line in list_of_line:

        ## Score and Expect
        if "Score" in line:
            line = line[:-1]  # remove "\n"
            S_line = string.split(line, " = ") 
            Expect = S_line[-1]  ## ***** Expect
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

            identities = S_line[-2]  ## ***** identities

            S_line2 = string.split(identities, "/")
            hits = string.atof(S_line2[0])                 ## ***** hits
            length_matched = string.atof(S_line2[1])       ## ***** length_matched
            abs_nb_differences = length_matched - hits     ## ***** abs_nb_differences

            identity_percent =   hits/length_matched * 100  ## ***** identity_percent

            divergence_percent = abs_nb_differences/length_matched*100  ## ***** divergence_percent

            ## treat gap line if any
            if g ==1:    # means there are gaps
                S_line3 = string.split(gaps_line, " ")
                gaps_part  = S_line3[-2]
                S_line4 = string.split(gaps_part, "/")
                gaps_number = string.atoi(S_line4[0])                                       ## ***** gaps_number
                
                real_differences = abs_nb_differences - gaps_number                         ## ***** real_differences
                real_divergence_percent = (real_differences/length_matched)*100             ## ***** real_divergence_percent
            else:
                gaps_number = 0               
                real_differences = 0
                real_divergence_percent = divergence_percent

        ## Frame
        elif "Frame" in line:
            line = line[:-1]  # remove "\n"
            S_line = string.split(line, " = ")
            frame = S_line[1]

    list_informations=[length_matched, Expect, Score, identities, hits, identity_percent, divergence_percent,gaps_number, real_divergence_percent, frame, length_matched]

    return(list_informations)

# Used in S06, S09, S11
def get_pairs(fasta_file_path):
    F2 = open(fasta_file_path, "r")
    list_pairwises = []
    while 1:
        next2 = F2.readline()
        if not next2:
            break
        if next2[0] == ">":
            fasta_name_query = next2[:-1]
            next3 = F2.readline()
            fasta_seq_query = next3[:-1]
            next3 = F2.readline()    ## jump one empty line (if any after the sequence)
            fasta_name_match = next3[:-1]
            next3 = F2.readline()
            fasta_seq_match = next3[:-1]
            pairwise = [fasta_name_query,fasta_seq_query,fasta_name_match,fasta_seq_match]
            
            ## ADD pairwise with condition
            list_pairwises.append(pairwise)
    F2.close()
    
    return(list_pairwises)

def extract_length(length_string):   # format length string = 57...902
    l3 = string.split(length_string, "...")
    n1 = string.atoi(l3[0])
    n2 = string.atoi(l3[1])
    length = n2-n1
    return(length)

def filter_redondancy(list_paireu, MIN_LENGTH):

    bash1 = {}
    list_pairout = []
    
    for pair in list_paireu:
         query_name = pair[0]
         query_seq = pair[1]
         match_name = pair[2]
         match_seq = pair[3]

         l1 = string.split(query_name, "||")
         short_query_name = l1[0][1:]
         length_matched =  extract_length(l1[1])
         l2 = string.split(match_name, "||")
         short_match_name = l2[0][1:]
         binom = "%s_%s" %(short_query_name, short_match_name)
         
         if binom not in bash1.keys():
             bash1[binom] = [query_name, query_seq, match_name, match_seq, length_matched]
         else:
             old_length = bash1[binom][-1]
             if length_matched > old_length:
                 bash1[binom] = [query_name, query_seq, match_name, match_seq, length_matched]

    
    for bino in bash1.keys():
        length = bash1[bino][-1]
        if length > MIN_LENGTH:
            list_pairout.append(bash1[bino])

    return(list_pairout)