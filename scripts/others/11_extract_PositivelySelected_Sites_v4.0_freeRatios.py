#!/usr/bin/python

## Author: Eric FONTANILLAS
## Date: 10.02.11
## Object: Extract likelihood of a codeml output runned in "BRANCH-MODEL"




############################
### Def 1: Extract data1 ###
############################
### For branch model (NULL MODEL ... no branch partition)
def extract_data_NULL(fileu):
    #print fileu
    
    while 1:
        next = fileu.readline()
        if not next:
            break


        if "lnL" in next:
            lnL = next[:-1]

            ff1 = string.split(lnL, ":")
            ff2 = ff1[2]
            ff3 = string.split(ff2, " ")
            ff4 = ff3[-1]
            np = ff4[:-1]
            np = string.atoi(np)

            ff5 = ff1[3]
            ff6 = string.split(ff5, " ")

            for element in ff6:
                if element == "":
                    index = ff6.index(element)
                    del ff6[index]

            lnL = float(ff6[0])

        if "omega (dN/dS) = " in next:
            lulu = string.split(next[:-1], " ")
            #print lulu
            w = float(lulu[-1])

        
    fileu.close()

    return(lnL, np, w)


#####################################################



############################
### Def 2: Extract data1 ### 
############################
### For branch model (2 categories of w)
def extract_data_ALTERNATE(fileu):
    borne1 = 0
    BASH={}
    
    while 1:
        next = fileu.readline()
        if not next:
            break


        if "lnL" in next:
            lnL = next[:-1]

            ## Get np
            ff1 = string.split(lnL, ":")
            ff2 = ff1[2]
            ff3 = string.split(ff2, " ")
            ff4 = ff3[-1]
            np = ff4[:-1]
            np = string.atoi(np)

            ## Get lnL
            ff5 = ff1[3]
            ff6 = string.split(ff5, " ")
            for element in ff6:
                if element == "":
                    index = ff6.index(element)
                    del ff6[index]
            lnL = float(ff6[0])
            
            ## Get list of branches (e.g. ["7..8", "7..1", ...]
            next = fileu.readline()    ## next line
            ff7 = string.split(next, " ")
            ff8 = []
            for element in ff7:
                if element != '' and element != '\n':
                    ff8.append(element)
            list_branches = ff8
            
            
        if "w (dN/dS) for branches" in next:
            borne1 = 1
            

        if borne1 == 1:
            for branch in list_branches:
                if branch in next:
                    #print next
                    ff10 = string.split(next[:-1], " ")
                    ff11 = []
                    for element in ff10:
                        if element != '' and element != '\n':
                            ff11.append(element)
                    #print ff11

                    key = ff11[0]

                    t = float(ff11[1])
                    N = float(ff11[2])
                    S = float(ff11[3])
                    w = float(ff11[4])
                    dN = float(ff11[5])
                    dS = float(ff11[6])
                    NdN = float(ff11[7])
                    SdS = float(ff11[8])
                    
                    values = [t,N,S,w,dN,dS,NdN,SdS]
                                      
                    BASH[key] = values

    fileu.close()
    return(lnL,np,list_branches, BASH)
#####################################################

###################
### RUN RUN RUN ###
###################

import string, os

## 1 ## Open input/output
path1 = "./TREE_FreeRatios"
path_tmp = "tmp"
L1 = os.listdir(path1)

os.mkdir("12_2_Branches")
branch_folder = "12_2_Branches"

## 2 ## Get list of sixtet names (pairs of NULL vs ALTERNATE sequences for each sixtet names)
listeu = []
for folder in L1:
    if "NULL" in folder:
        name = folder[:-5]
        #print name
        listeu.append(name)


## 3 ## Write header for general table ==>  open 1rst folder (i.e. locus) in order to get the HEADERS  
OUT1 = open("12_1_out.csv", "w")
OUT1.write("Name;lnL_NULL;lnL_ALTERNATE;df;NULL_w;")
sixtet_name1 = listeu[0]

#path_Null = "%s/%s_NULL" %(path1, sixtet_name1)
path_Alternate = "%s/%s_ALTERNATE" %(path1, sixtet_name1)
#IN_N = open("%s/sixtet.out" %path_Null, "r")
IN_A = open("%s/sixtet.out" %path_Alternate, "r") # open 1rst folder (i.e. locus) in order to get the HEADERS 

A_lnL, A_np, list_branches,BASH = extract_data_ALTERNATE(IN_A)  ### DEF 2 ###
for branch in list_branches:
    t = "t_" + branch
    N = "N_" + branch
    S = "S_" + branch
    w = "w_" + branch
    dN = "dN_" + branch
    dS = "dS_" + branch
    NdN = "NdN_" + branch
    SdS = "SdS_" + branch
    NdN_std = "NdN/(NdN+SdS)_" + branch
    SdS_std = "SdS/(NdN+SdS)_" + branch
    OUT1.write("%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;" %(w,t,dN,dS,NdN,SdS,NdN_std,SdS_std,N,S))
OUT1.write("\n")

## 4 ## Write header for INDIVIDUAL tables (1 table per branch)
for branch in list_branches:
    ## Open individual output
    OUT_branch = open("%s/%s.csv" %(branch_folder,branch), "w")
    ## Print header in individual output
    OUT_branch.write("Locus_name;w;t;dN;dS;NdN;SdS;NdN/(NdN+SdS);SdS/(NdN+SdS);N;S\n")
    OUT_branch.close()
    
## 5 ## Explore results for each locus
for sixtet_name in listeu:
    OUT1.write("%s;" %sixtet_name)
    path_Null = "%s/%s_NULL" %(path1, sixtet_name)
    path_Alternate = "%s/%s_ALTERNATE" %(path1, sixtet_name)
    IN_N = open("%s/sixtet.out" %path_Null, "r")
    IN_A = open("%s/sixtet.out" %path_Alternate, "r")
    N_lnL, N_np, N_w = extract_data_NULL(IN_N)   ### DEF 1 ###
    A_lnL, A_np, list_branches,BASH = extract_data_ALTERNATE(IN_A)  ### DEF 2 ###

    ## Compute df
    df = A_np - N_np

    ## Print output    
    OUT1.write("%.5f;%.5f;%d;%.5f;" %(N_lnL,A_lnL,df,N_w))

    for branch in list_branches:
        values = BASH[branch]

        ## Re-Open individual output with "ADD" option!!!
        OUT_branch = open("%s/%s.csv" %(branch_folder,branch), "a")

        ## Get the values per branch
        values = BASH[branch]

        t = values[0]
        N  = values[1]
        S = values[2]
        w = values[3]
        dN = values[4]
        dS = values[5]
        NdN = values[6]
        SdS = values[7]
        
        SUM = NdN + SdS

        if t > 0 and SUM > 0:
            NdN_std = float(NdN)/float(SUM)
            SdS_std = float(SdS)/float(SUM)
        else:
            NdN_std = 0
            SdS_std = 0


        OUT1.write("%.5f;%.5f;%.5f;%.5f;%.5f;%.5f;%.5f;%.5f;%.5f;%.5f;" %(w,t,dN,dS,NdN,SdS,NdN_std,SdS_std,N,S))
        OUT_branch.write("%s;%.5f;%.5f;%.5f;%.5f;%.5f;%.5f;%.5f;%.5f;%.5f;%.5f\n" %(sixtet_name,w,t,dN,dS,NdN,SdS,NdN_std,SdS_std,N,S))
        OUT_branch.close()

    OUT1.write("\n")
        
    
    IN_N.close()
    IN_A.close()

OUT1.close()
