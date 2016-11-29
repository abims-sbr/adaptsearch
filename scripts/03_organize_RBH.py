#!/usr/local/public/python-2.6.5/bin/python2.6

#/usr/bin/python

## AUTHOR: Eric Fontanillas

## LAST VERSION: 23.05.2011

## DESCRIPTION: 


import sys, string, os, itertools

#PATH_PATRON_BLASTN = "02_PATRON_RBH_BLASTN"
PATH_PATRON_TBLASTX = "02_PATRON_RBH_TBLASTX"

PATH_IN_TMP = "../tmp/01_formated_INPUT"
PATH_RBH_TMP = "../tmp/02_ParwiseRBH_tblastx"

L0 = os.listdir(PATH_RBH_TMP)
if L0 != []:
    os.system("rm -fr %s/*" %PATH_RBH_TMP)

L1 = os.listdir(PATH_IN_TMP)

MULTI_QSUB = open("04_multiQSUB.sh", "w")



##############################################
list_pairwise = itertools.combinations(L1, 2)
##############################################

for pairwise in list_pairwise:
    print pairwise
    DIR1 = pairwise[0]
    DIR2 = pairwise[1]

    path_DIR1 = "%s/%s" %(PATH_IN_TMP, DIR1)
    path_DIR2 = "%s/%s" %(PATH_IN_TMP, DIR2)

    
    
    S1 = string.split(DIR1, "_")
    shortDIR1 = S1[0]

    S2 = string.split(DIR2, "_")
    shortDIR2 = S2[0]

    RBH_folder = shortDIR1 + "_" + shortDIR2
    path_RBH_folder = PATH_RBH_TMP + "/" + RBH_folder

    os.mkdir(path_RBH_folder)

    os.system("cp -fr %s/0* %s/" %(PATH_PATRON_TBLASTX, path_RBH_folder))
    os.system("cp -fr %s/1* %s/" %(PATH_PATRON_TBLASTX, path_RBH_folder))
    os.system("cp -fr %s/2* %s/" %(PATH_PATRON_TBLASTX, path_RBH_folder))
    os.system("cp -fr %s/Ortho.qsub %s/XOrtho_%s_%s.qsub" %(PATH_PATRON_TBLASTX, path_RBH_folder, shortDIR1, shortDIR2))
    
    os.system("cp -fr %s %s/"  %(path_DIR1, path_RBH_folder))
    os.system("cp -fr %s %s/"  %(path_DIR2, path_RBH_folder))
    
    pipeline_script = open("%s/XXX_pipeline.sh" %path_RBH_folder, "w")
    pipeline_patron = open("%s/XXX_patronPipeline.sh" %PATH_PATRON_TBLASTX, "r")

    F1 = pipeline_patron.read()
    F1 = string.replace(F1, "XX1", DIR1)
    F1 = string.replace(F1, "XX2", DIR2)
    
    pipeline_script.write(F1)

    pipeline_patron.close()
    pipeline_script.close()

    os.system("chmod +x %s/*" %path_RBH_folder)

    os.mkdir("%s/tmp" %path_RBH_folder)


    MULTI_QSUB.write("cd %s\n" %path_RBH_folder)
    MULTI_QSUB.write("qsub -q long.q XOrtho_%s_%s.qsub\n" %(shortDIR1, shortDIR2))
    MULTI_QSUB.write("cd ../../../script/\n")
    

MULTI_QSUB.close()



os.system("chmod +x 04_multiQSUB.sh")
