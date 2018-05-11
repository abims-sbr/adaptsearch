#!/usr/bin/env python
#coding: utf-8
#Author : Eric Fontanillas (2010) - Victor Mataigne (2018)

import os
#import pandas as pd

def tableu(fileu):
    F1 = open(fileu, "r")
    basheu = {}
    
    ## jump headers (and get them, but not usefull here)
    HEADERS = F1.readline()
    S = HEADERS[:-1].split(",")
    list_HEADERS = S[1:]
    #print list_HEADERS
    
    # Initialize counts
    Pp_greater = 0
    Pp_lower = 0

    Pg_greater = 0
    Pg_lower = 0
    
    Ap_greater = 0
    Ap_lower = 0

    Ps_greater = 0
    Ps_lower = 0

    while 1:
        nextline = F1.readline()
        if not nextline :
            break
        ## Get values
        S = nextline[:-1].split(",")
        locus_name = S[0]
        Ac = float(S[1])
        Ap = float(S[2])
        Pf = float(S[3])
        Pg = float(S[4])
        Pp = float(S[5])
        Ps = float(S[6])
        
        ## Get max HOTS and min HOTS
        MAX_ApPs = max(Ap,Ps)
        MIN_ApPs = min(Ap,Ps)
        
        ## Get max COLDS and min COLDS
        MAX_PpPg = max(Pp,Pg)
        MIN_PpPg = min(Pp,Pg)

        ## TESTS : Pp vs HOTS
        if Pp > MAX_ApPs:
            Pp_greater = Pp_greater + 1
        elif Pp < MIN_ApPs:
            Pp_lower = Pp_lower + 1

        ## TESTS : Pg vs HOTS
        if Pg > MAX_ApPs:
            Pg_greater = Pg_greater + 1
        elif Pg < MIN_ApPs:
            Pg_lower = Pg_lower + 1

        ## TESTS : Ap vs COLDS
        if Ap > MAX_PpPg:
            Ap_greater = Ap_greater + 1
        elif Ap < MIN_PpPg:
            Ap_lower = Ap_lower + 1

        ## TESTS : Ps vs COLDS
        if Ps > MAX_PpPg:
            Ps_greater = Ps_greater + 1
        elif Ps < MIN_PpPg:
            Ps_lower = Ps_lower + 1

    return(Pp_greater, Pp_lower, Pg_greater, Pg_lower, Ap_greater, Ap_lower, Ps_greater, Ps_lower)
#####################

#def tableu_on_aa():

#def tableu_on_variables():

def main():

    LAA =['K','R','A','F','I','L','M','V','W','N','Q','S','T','H','Y','C','D','E','P','G']
    LV = ['IVYWREL','EK','ERK','DNQTSHA','QH','ratio_ERK_DNQTSHA','ratio_EK_QH','FYMINK','GARP',
          'ratio_GARP_FYMINK','AVLIM','FYW','AVLIMFYW','STNQ','RHK','DE','RHKDE','APGC','AC',
          'VLIM','ratio_AC_VLIM','ratio_APGC_VLIM']

    input_path_aa = '02_tables_per_aa'
    list_in_aa = os.listdir(input_path_aa)
    input_path_var = '02_tables_per_variable'
    list_in_var = os.listdir(input_path_var)
    out_path_aa = '03_tables_counts_signTest_aa'
    os.mkdir(out_path_aa)
    out_path_var = '03_tables_counts_signTest_variables'
    os.mkdir(out_path_var)

    for variable in LAA:
        F_IN = "%s/%s" %(input_path_aa, variable)
        F_OUT = open("%s/%s.csv" %(out_path_aa, variable), "w")

        Pp_greater, Pp_lower, Pg_greater, Pg_lower, Ap_greater, Ap_lower, Ps_greater, Ps_lower = tableu(F_IN)    ### DEF 1 ###       

        Pp_diff = Pp_greater - Pp_lower
        Pg_diff = Pg_greater - Pg_lower
        Ap_diff = Ap_greater - Ap_lower
        Ps_diff = Ps_greater - Ps_lower        

        Pp_total = Pp_greater + Pp_lower
        Pg_total = Pg_greater + Pg_lower
        Ap_total = Ap_greater + Ap_lower
        Ps_total = Ps_greater + Ps_lower

        headers = ",Pp_vs_HOTS,Pg_vs_HOTS,Ap_vs_COLDS,Ps_vs_COLDS\n"
        line_greater = "Greater,%d,%d,%d,%d\n" %(Pp_greater,Pg_greater,Ap_greater,Ps_greater)
        line_lower = "Lower,%d,%d,%d,%d\n" %(Pp_lower,Pg_lower,Ap_lower,Ps_lower)
        line_diff = "Difference,%d,%d,%d,%d\n" %(Pp_diff,Pg_diff,Ap_diff,Ps_diff)
        line_total = "Trial number,%d,%d,%d,%d\n" %(Pp_total,Pg_total,Ap_total,Ps_total)
        
        F_OUT.write(headers)
        F_OUT.write(line_greater)
        F_OUT.write(line_lower)
        F_OUT.write(line_diff)
        F_OUT.write(line_total)
        F_OUT.close()

    for variable in LV:        
        F_IN = "%s/%s" %(input_path_var, variable)
        F_OUT = open("%s/%s.csv" %(out_path_var, variable), "w")

        Pp_greater, Pp_lower, Pg_greater, Pg_lower, Ap_greater, Ap_lower, Ps_greater, Ps_lower = tableu(F_IN)    ### DEF 1 ###

        Pp_diff = Pp_greater - Pp_lower
        Pg_diff = Pg_greater - Pg_lower
        Ap_diff = Ap_greater - Ap_lower
        Ps_diff = Ps_greater - Ps_lower        

        Pp_total = Pp_greater + Pp_lower
        Pg_total = Pg_greater + Pg_lower
        Ap_total = Ap_greater + Ap_lower
        Ps_total = Ps_greater + Ps_lower

        headers = ",Pp_vs_HOTS,Pg_vs_HOTS,Ap_vs_COLDS,Ps_vs_COLDS\n"
        line_greater = "Greater,%d,%d,%d,%d\n" %(Pp_greater,Pg_greater,Ap_greater,Ps_greater)
        line_lower = "Lower,%d,%d,%d,%d\n" %(Pp_lower,Pg_lower,Ap_lower,Ps_lower)
        line_diff = "Difference,%d,%d,%d,%d\n" %(Pp_diff,Pg_diff,Ap_diff,Ps_diff)
        line_total = "Trial number,%d,%d,%d,%d\n" %(Pp_total,Pg_total,Ap_total,Ps_total)
        
        F_OUT.write(headers)
        F_OUT.write(line_greater)
        F_OUT.write(line_lower)
        F_OUT.write(line_diff)
        F_OUT.write(line_total)
        F_OUT.close()

    os.mkdir('04_outputs_aa')
    os.mkdir('04_outputs_variables')
    os.system('Rscript new_sign_test_binomial.R --indir %s --outdir 04_outputs_aa' %out_path_aa)
    os.system('Rscript new_sign_test_binomial.R --indir %s --outdir 04_outputs_variables' %out_path_var)

if __name__ == '__main__':
    main()