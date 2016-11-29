#!/usr/bin/python


MINIMUM_LENGTH = 1


############################
##### DEF1 : Get Pairs #####
############################
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
            #next3 = F2.readline()
            fasta_name_match = next3[:-1]
            next3 = F2.readline()
            fasta_seq_match = next3[:-1]
            pairwise = [fasta_name_query,fasta_seq_query,fasta_name_match,fasta_seq_match]
            
            ## ADD pairwise with condition
            list_pairwises.append(pairwise)
    F2.close()
    #print list_pairwises
    return(list_pairwises)
##############################################

#################################
##### DEF2 : Extract length #####
#################################
def extract_length(length_string):   # format length string = 57...902
    l3 = string.split(length_string, "...")
    n1 = string.atoi(l3[0])
    n2 = string.atoi(l3[1])
    length = n2-n1
    return(length)

##############################################

####################################
##### DEF3 : Remove Redondancy #####
####################################
def filter_redondancy(list_paireu, MIN_LENGTH):

    bash1 = {}
    list_pairout = []
    
    for pair in list_paireu:
         query_name = pair[0]
         #print query_name
         query_seq = pair[1]
         match_name = pair[2]
         match_seq = pair[3]

         l1 = string.split(query_name, "||")
         short_query_name = l1[0][1:]
         length_matched =  extract_length(l1[1])          ### DEF2 ###
         l2 = string.split(match_name, "||")
         short_match_name = l2[0][1:]
         binom = "%s_%s" %(short_query_name, short_match_name)

         #print binom
         
         if binom not in bash1.keys():
             bash1[binom] = [query_name, query_seq, match_name, match_seq, length_matched]
         else:
             old_length = bash1[binom][-1]
             if length_matched > old_length:
                 bash1[binom] = [query_name, query_seq, match_name, match_seq, length_matched]
             #else:
             #    print "shorter!"

    
    for bino in bash1.keys():
        print bino
        length = bash1[bino][-1]
        #print length
        if length > MIN_LENGTH:
            print "\t %d bp" %(length)
            list_pairout.append(bash1[bino])
        else:
            print "\t%s have been removed due to length (%d bp)" %(bino, length)


    #print bash1.keys()

    return(list_pairout)

##############################################


################################
###### DEF3: filter lengh ######
################################
# def filter_length(list_pairs, MIN_LENGTH):
#     list_pairs2 = []
#     for pair in list_pairs:
#         #print pair[0]
#         seq = pair[1]
#         #print seq
#         ln = len(seq)
#         #print ln
#         if ln > MIN_LENGTH:
#             list_pairs2.append(pair)
#     print list_pairs2
#     return(list_pairs2)
##############################################

#######################
##### RUN RUN RUN #####
#######################
import string, os, time, re, sys

WORK_DIR = sys.argv[1]

F_IN = "%s/17_ReciprocalHits.fasta" %WORK_DIR
F_OUT = "%s/19_ReciprocalBestHits.fasta"  %(WORK_DIR)
File_OUT = open(F_OUT, "w")
F_OUT2 = "%s/19_onlyQuery_RBH_ungapped.fasta"  %(WORK_DIR)
File_OUT2 = open(F_OUT2, "w")
F_OUT3 = "%s/19_onlyQuery_RBH_ungapped_shortName.fasta"  %(WORK_DIR)
File_OUT3 = open(F_OUT3, "w")

list_pairwises = get_pairs(F_IN)                                                           ### DEF1 ###
#for paireu in list_pairwises:
#    print paireu[0]

list_pairwises_filtered1 = filter_redondancy(list_pairwises, MINIMUM_LENGTH)                               ### DEF3 ###



i = 0
for pair in list_pairwises_filtered1:
     i = i+1

     ## Write pairwise alignment
     File_OUT.write("%s\n" %pair[0])
     File_OUT.write("%s\n" %pair[1])
     File_OUT.write("%s\n" %pair[2])
     File_OUT.write("%s\n" %pair[3])
     
     seq = pair[1]
     ungapped_seq = string.replace(seq, "-", "")

     name = pair[0]
     S1 = string.split(name, "||")
     short_name = S1[0]
     
     File_OUT2.write("%s\n" %name)
     File_OUT2.write("%s\n" %ungapped_seq)
     
     File_OUT3.write("%s\n" %short_name)
     File_OUT3.write("%s\n" %ungapped_seq)     
     
File_OUT.close()
File_OUT2.close()
File_OUT3.close()
##############################################
