#!/usr/bin/python

## Author: Eric FONTANILLAS
## Date: 21.12.10
## Object: Test for compositional bias in genome and proteome as marker of thermal adaptation (comparison between 2 "hot" species: Ap and Ps and one "cold" species: Pg)

#############
### DEF 0 ###
#############
def simplify_fasta_name(fasta_name):
    L=["Ap", "Ac", "Ps", "Pf", "Pg", "Pp"]

    for abbreviation in L:
        if abbreviation in fasta_name:
            new_fasta_name = abbreviation

    return(new_fasta_name)
##########################################

###########
## DEF1 ##
###########
## Generates bash, with key = fasta name; value = sequence (WITH GAP, IF ANY, REMOVED IN THIS FUNCTION)

def dico(fasta_file):

    count_fastaName=0
    F1 = open(fasta_file, "r")
    
    bash1 = {}
    while 1:
        nextline = F1.readline()
        #print nextline
        if not nextline :
            break
        
        if nextline[0] == ">":
            count_fastaName = count_fastaName + 1
            fasta_name = nextline[1:-1]
            nextline = F1.readline()
            sequence = nextline[:-1]
            
            if fasta_name not in bash1.keys():
                fasta_name = simplify_fasta_name(fasta_name)  ### DEF 0 ###
                bash1[fasta_name] = sequence
            else:
                print fasta_name

    # Find alignment length
    kk = bash1.keys()
    key0 = kk[0]
    seq0 = bash1[key0]
    ln_seq = len(seq0)

    F1.close()
    
    return(bash1)
#####################################



##################
###### DEF2 ######
##################
def base_composition(seq):
      count_A=string.count(seq, "A")
      count_T=string.count(seq, "T")
      count_C=string.count(seq, "C")
      count_G=string.count(seq, "G")


      CG = count_C+count_G
      AT = count_T+count_A
      
      AG = count_A+count_G
      TC = count_T+count_C
      
      ## 1 ## Search for compositional bias in genome as marker of thermal adaptation: CG vs AT
      ratio_CG_AT=float(CG)/float(AT)
      
      ## 2 ## Search for compositional bias in genome as marker of thermal adaptation: AG vs TC
      ratio_purine_pyrimidine=float(AG)/float(TC)

      ## 3 ## Nucleotide proportion
      ln = len(seq)
      prop_A = float(count_A)/float(ln)
      prop_T = float(count_T)/float(ln)
      prop_C = float(count_C)/float(ln)
      prop_G = float(count_G)/float(ln)
      

      return(ratio_CG_AT, ratio_purine_pyrimidine, prop_A, prop_T, prop_C, prop_G)
##############################################


##################
###### DEF3 ######
##################
def aa_composition1(seq):
    
    ## 1 ## count occurence of AA
    count_K=string.count(seq,"K")
    count_R=string.count(seq,"R")
    count_A=string.count(seq,"A")
    count_F=string.count(seq,"F")
    count_I=string.count(seq,"I")
    count_L=string.count(seq,"L")
    count_M=string.count(seq,"M")
    count_V=string.count(seq,"V")
    count_W=string.count(seq,"W")
    count_N=string.count(seq,"N")
    count_Q=string.count(seq,"Q")
    count_S=string.count(seq,"S")
    count_T=string.count(seq,"T")
    count_H=string.count(seq,"H")
    count_Y=string.count(seq,"Y")
    count_C=string.count(seq,"C")
    count_D=string.count(seq,"D")
    count_E=string.count(seq,"E")
    count_P=string.count(seq,"P")
    count_G=string.count(seq,"G")

    ## 2 ## compute relative proportion

    ln = len(seq)
    
    prop_K=float(count_K)/float(ln)
    prop_R=float(count_R)/float(ln)
    prop_A=float(count_A)/float(ln)
    prop_F=float(count_F)/float(ln)
    prop_I=float(count_I)/float(ln)
    prop_L=float(count_L)/float(ln)
    prop_M=float(count_M)/float(ln)
    prop_V=float(count_V)/float(ln)
    prop_W=float(count_W)/float(ln)
    prop_N=float(count_N)/float(ln)
    prop_Q=float(count_Q)/float(ln)
    prop_S=float(count_S)/float(ln)
    prop_T=float(count_T)/float(ln)
    prop_H=float(count_H)/float(ln)
    prop_Y=float(count_Y)/float(ln)
    prop_C=float(count_C)/float(ln)
    prop_D=float(count_D)/float(ln)
    prop_E=float(count_E)/float(ln)
    prop_P=float(count_P)/float(ln)
    prop_G=float(count_G)/float(ln)


        
    return(prop_K,prop_R,prop_A,prop_F,prop_I,prop_L,prop_M,prop_V,prop_W,prop_N,prop_Q,prop_S,prop_T,prop_H,prop_Y,prop_C,prop_D,prop_E,prop_P,prop_G)
      
##################
###### DEF4 ######
##################
def aa_composition2(seq):
    
    ## 1 ## count occurence of AA
    count_K=string.count(seq,"K")
    count_R=string.count(seq,"R")
    count_A=string.count(seq,"A")
    count_F=string.count(seq,"F")
    count_I=string.count(seq,"I")
    count_L=string.count(seq,"L")
    count_M=string.count(seq,"M")
    count_V=string.count(seq,"V")
    count_W=string.count(seq,"W")
    count_N=string.count(seq,"N")
    count_Q=string.count(seq,"Q")
    count_S=string.count(seq,"S")
    count_T=string.count(seq,"T")
    count_H=string.count(seq,"H")
    count_Y=string.count(seq,"Y")
    count_C=string.count(seq,"C")
    count_D=string.count(seq,"D")
    count_E=string.count(seq,"E")
    count_P=string.count(seq,"P")
    count_G=string.count(seq,"G")

    ## 2 ## compute seq length

    ln = len(seq)

    ## 3 ## Famous Hyperthermophile Prokaryotes criterias
    
    # 3.1. IVYWREL estimator  => positivelly correlated with otpimal growth
    count_IVYWREL = count_I+count_V+count_Y+count_W+count_R+count_E+count_L
    prop_IVYWREL = float(count_IVYWREL)/float(ln)
    
    # 3.2. ERK estimator (i.e. ERK vs DNQTSHA) => positivelly correlated with optimal growth temperature
    # ERK alone
    count_ERK = count_E + count_R + count_K
    prop_ERK = float(count_ERK)/float(ln)
    # DNQTSHA alone
    count_DNQTSHA = count_D + count_N + count_Q +  count_T + count_S + count_H + count_A
    prop_DNQTSHA=float(count_DNQTSHA)/float(ln)
    # ERK vs DNQTSHA
    if count_DNQTSHA != 0:
        ratio_ERK_vs_DNQTSHA = float(count_ERK)/float(count_DNQTSHA)
    else:
        ratio_ERK_vs_DNQTSHA=-1
    # EK/QH estimator
    count_EK = count_E+count_K
    count_QH = count_Q+count_H

    prop_EK = float(count_EK)/float(ln)
    prop_QH = float(count_QH)/float(ln)

    if count_QH != 0:
        ratio_EK_vs_QH = float(count_EK)/float(count_QH)
    else:
        ratio_EK_vs_QH=-1     ## "-1" will indicate the impossibility to compute the ratio (coz the numerator)

    ## 4 ## Mutationnal bias hypothesis => AT rich: favor FYMINK   // GC rich: favor GARP
    ## The mutational bias model predict a linear relationship between GARP vs FYMINK    ==> so if outliers to that, it means that the excess of GARP or FYMINK are not explained by the mutationnal bias model but by other thing ... selection!!???
    count_FYMINK=count_F+count_Y+count_M+count_I+count_N+count_K
    prop_FYMINK = float(count_FYMINK)/float(ln)

    count_GARP=count_G+count_A+count_R+count_P
    prop_GARP=float(count_GARP)/float(ln)

    ## 5 ## Hydophobicity hypothesis [should INCREASE with thermal adaptation]
    ## 5.1. AL
    count_AVLIMFYW = count_A+count_V+count_L+count_I+count_F+count_Y+count_W+count_M
    prop_AVLIMFYW=float(count_AVLIMFYW)/float(ln)
    ## 5.2. Only non-aromatic
    count_AVLIM = count_A+count_V+count_L+count_I+count_M
    prop_AVLIM=float(count_AVLIM)/float(ln)
    ## 5.3. Only aromatic (have they higher residus volume?? in such case opposite hypothesis based on residu volume, predict DECREASE for these aa in composition) 
    count_FYW = count_F+count_Y+count_W
    prop_FYW=float(count_FYW)/float(ln)

    ## 6 ## Charged hypothesis  => positivelly correlated with optimal growth temperature
    # All charged
    count_RHKDE = count_R + count_H +count_K + count_D + count_E
    prop_RHKDE = float(count_RHKDE)/float(ln)
    # Only positive
    count_RHK = count_R + count_H +count_K
    prop_RHK = float(count_RHK)/float(ln)
    # Only negative
    count_DE = count_D + count_E
    prop_DE = float(count_DE)/float(ln)
    
    ## 7 ## Neutral polar hypothesis  [should DECREASE with thermal adaptation]
    count_STNQ = count_S+count_T+count_N+count_Q
    prop_STNQ=float(count_STNQ)/float(ln)


    ## 9 ## Didier's Theory ==> advantage in mutation for smaller residues in the most common mutation (excluding charged and aromatique), for thermophyle.
    ## 9.1 ## Didier's criteria 1 = SMALL / BIG
    count_APGC = count_A+count_G+count_P+count_C
    prop_APGC=float(count_APGC)/float(ln)
    count_VLIM = count_V+count_L+count_I+count_M
    prop_VLIM=float(count_VLIM)/float(ln)
    if count_VLIM != 0:
        ratio_APGC_vs_VLIM = float(count_APGC)/float(count_VLIM)
    else:
        ratio_APGC_vs_VLIM=-1     ## "-1" will indicate the impossibility to compute the ratio (coz the numerator)

    ## 9.2 ## Didier's criteria 2 = VERY SMALL / BIG
    count_AC = count_A+count_C
    prop_AC=float(count_AC)/float(ln)

    #count_VLIM = count_V+count_L+count_I+count_M
    if count_VLIM != 0:
        ratio_AC_vs_VLIM = float(count_AC)/float(count_VLIM)
    else:
        ratio_AC_vs_VLIM=-1     ## "-1" will indicate the impossibility to compute the ratio (coz the numerator)

    return(count_IVYWREL,prop_IVYWREL,count_ERK,prop_ERK,count_DNQTSHA,prop_DNQTSHA,ratio_ERK_vs_DNQTSHA,count_EK,prop_EK,count_QH,prop_QH,ratio_EK_vs_QH,count_FYMINK,prop_FYMINK,count_GARP,prop_GARP,count_AVLIMFYW, prop_AVLIMFYW,count_AVLIM,prop_AVLIM,count_FYW,prop_FYW,count_STNQ, prop_STNQ, count_VLIM,prop_VLIM, count_APGC,prop_APGC, count_AC,prop_AC, ratio_APGC_vs_VLIM, ratio_AC_vs_VLIM, count_RHKDE,prop_RHKDE,count_RHK,prop_RHK,count_DE,prop_DE)
#####################


##################
###### DEF5 ######
##################
def aa_properties(fileIN_aaProperties):
    next = fileIN_aaProperties.readline()    ## JUMP HEADERS

    bash_aa_properties={}

    while 1:
        next = fileIN_aaProperties.readline()
        if not next:
            break

        S1 = string.split(next, ",")

        aa_name = S1[1]
        S2 = string.split(aa_name, "/")
        aa_code = S2[1][:-1]

        frequencies = S1[2][:-1]
        Residue_Weight = S1[5]
        Residue_Volume = S1[6]
        Partial_specific_volume = S1[7]
        Hydration = S1[8]

        bash_aa_properties[aa_code] = [frequencies,Residue_Weight,Residue_Volume,Partial_specific_volume,Hydration]

    return(bash_aa_properties)
    
    
##################
###### DEF6 ######
##################
def sequence_properties_from_aa_properties(seq, bash_properties):
    
    ## 1 ## count occurence of AA
    count_K=string.count(seq,"K")
    count_R=string.count(seq,"R")
    count_A=string.count(seq,"A")
    count_F=string.count(seq,"F")
    count_I=string.count(seq,"I")
    count_L=string.count(seq,"L")
    count_M=string.count(seq,"M")
    count_V=string.count(seq,"V")
    count_W=string.count(seq,"W")
    count_N=string.count(seq,"N")
    count_Q=string.count(seq,"Q")
    count_S=string.count(seq,"S")
    count_T=string.count(seq,"T")
    count_H=string.count(seq,"H")
    count_Y=string.count(seq,"Y")
    count_C=string.count(seq,"C")
    count_D=string.count(seq,"D")
    count_E=string.count(seq,"E")
    count_P=string.count(seq,"P")
    count_G=string.count(seq,"G")


    ## 2 ## Compute properties 1: Residue Weight (Mr) (UNIT:Daltons):

    Total_Residue_Weight = count_K*float(bash_properties["K"][1]) + count_R*float(bash_properties["R"][1]) + count_A*float(bash_properties["A"][1]) + count_F*float(bash_properties["F"][1]) + count_I*float(bash_properties["I"][1]) + count_L*float(bash_properties["L"][1]) + count_M*float(bash_properties["M"][1]) + count_V*float(bash_properties["V"][1]) + count_W*float(bash_properties["W"][1]) + count_N*float(bash_properties["N"][1]) + count_Q*float(bash_properties["Q"][1]) + count_S*float(bash_properties["S"][1]) + count_T*float(bash_properties["T"][1]) + count_H*float(bash_properties["H"][1]) + count_Y*float(bash_properties["Y"][1]) + count_C*float(bash_properties["C"][1]) + count_D*float(bash_properties["D"][1]) + count_E*float(bash_properties["E"][1]) + count_P*float(bash_properties["P"][1]) + count_G*float(bash_properties["G"][1])
    Total_Residue_Volume = count_K*float(bash_properties["K"][2]) + count_R*float(bash_properties["R"][2]) + count_A*float(bash_properties["A"][2]) + count_F*float(bash_properties["F"][2]) + count_I*float(bash_properties["I"][2]) + count_L*float(bash_properties["L"][2]) + count_M*float(bash_properties["M"][2]) + count_V*float(bash_properties["V"][2]) + count_W*float(bash_properties["W"][2]) + count_N*float(bash_properties["N"][2]) + count_Q*float(bash_properties["Q"][2]) + count_S*float(bash_properties["S"][2]) + count_T*float(bash_properties["T"][2]) + count_H*float(bash_properties["H"][2]) + count_Y*float(bash_properties["Y"][2]) + count_C*float(bash_properties["C"][2]) + count_D*float(bash_properties["D"][2]) + count_E*float(bash_properties["E"][2]) + count_P*float(bash_properties["P"][2]) + count_G*float(bash_properties["G"][2])
    Total_Partial_specific_volume = count_K*float(bash_properties["K"][3]) + count_R*float(bash_properties["R"][3]) + count_A*float(bash_properties["A"][3]) + count_F*float(bash_properties["F"][3]) + count_I*float(bash_properties["I"][3]) + count_L*float(bash_properties["L"][3]) + count_M*float(bash_properties["M"][3]) + count_V*float(bash_properties["V"][3]) + count_W*float(bash_properties["W"][3]) + count_N*float(bash_properties["N"][3]) + count_Q*float(bash_properties["Q"][3]) + count_S*float(bash_properties["S"][3]) + count_T*float(bash_properties["T"][3]) + count_H*float(bash_properties["H"][3]) + count_Y*float(bash_properties["Y"][3]) + count_C*float(bash_properties["C"][3]) + count_D*float(bash_properties["D"][3]) + count_E*float(bash_properties["E"][3]) + count_P*float(bash_properties["P"][3]) + count_G*float(bash_properties["G"][3])
    Total_Hydration = count_K*float(bash_properties["K"][4]) + count_R*float(bash_properties["R"][4]) + count_A*float(bash_properties["A"][4]) + count_F*float(bash_properties["F"][4]) + count_I*float(bash_properties["I"][4]) + count_L*float(bash_properties["L"][4]) + count_M*float(bash_properties["M"][4]) + count_V*float(bash_properties["V"][4]) + count_W*float(bash_properties["W"][4]) + count_N*float(bash_properties["N"][4]) + count_Q*float(bash_properties["Q"][4]) + count_S*float(bash_properties["S"][4]) + count_T*float(bash_properties["T"][4]) + count_H*float(bash_properties["H"][4]) + count_Y*float(bash_properties["Y"][4]) + count_C*float(bash_properties["C"][4]) + count_D*float(bash_properties["D"][4]) + count_E*float(bash_properties["E"][4]) + count_P*float(bash_properties["P"][4]) + count_G*float(bash_properties["G"][4])
    
    
    return(Total_Residue_Weight,Total_Residue_Volume,Total_Partial_specific_volume,Total_Hydration)

########################################################



###################
### RUN RUN RUN ###
###################
import string, os


## 1 ## List taxa
LT = ["Ap","Ac","Pg","Pf","Ps", "Pp"]

## 2 ## PathIN
fileIN_properties = open("01_AminoAcid_Properties2.csv", "r")

# Path_IN_loci_NUC = "02_CDS_No_Missing_Data_nuc_CDS_withM"
# Lloci_NUC = os.listdir(Path_IN_loci_NUC)

Path_IN_loci_AA = "02_3_concatenation_PolymorphicSites"
#Path_IN_loci_AA = "02_CDS_No_Missing_Data_aa_CDS_withM"
Lloci_AA = os.listdir(Path_IN_loci_AA)

## 3 ## PathOUT

## 3.1 ## PROT composition
fileOUT_PROT_ALL=open("13_prot_compositions_All_AA.csv","w")
fileOUT_PROT_ALL.write("LOCUS,")
for taxa in LT:
    fileOUT_PROT_ALL.write("%s_prop_K,%s_prop_R,%s_prop_A,%s_prop_F,%s_prop_I,%s_prop_L,%s_prop_M,%s_prop_V,%s_prop_W,%s_prop_N,%s_prop_Q,%s_prop_S,%s_prop_T,%s_prop_H,%s_prop_Y,%s_prop_C,%s_prop_D,%s_prop_E,%s_prop_P,%s_prop_G," %(taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa))
fileOUT_PROT_ALL.write("\n")

## 3.2 ## PROT IVYWREL
fileOUT_IVYWREL=open("14_IVYWREL.csv","w")
fileOUT_IVYWREL.write("LOCUS,")
for taxa in LT:
    fileOUT_IVYWREL.write("%s_count_IVYWREL,%s_prop_IVYWREL," %(taxa,taxa))
fileOUT_IVYWREL.write("\n")


## 3.3 ## PROT ERK_DNQTSHA
fileOUT_ERK_DNQTSHA=open("15_ERK_DNQTSHA.csv","w")
fileOUT_ERK_DNQTSHA.write("LOCUS,")
for taxa in LT:
    fileOUT_ERK_DNQTSHA.write("%s_count_ERK,%s_prop_ERK,%s_count_DNQTSHA,%s_prop_DNQTSHA,%s_ratio_ERK_vs_DNQTSHA," %(taxa,taxa,taxa,taxa,taxa))
fileOUT_ERK_DNQTSHA.write("\n")

## 3.4 ## PROT EK_QH
fileOUT_EK_QH=open("16_EK_QH.csv","w")
fileOUT_EK_QH.write("LOCUS,")
for taxa in LT:
    fileOUT_EK_QH.write("%s_count_EK,%s_prop_EK,%s_count_QH,%s_prop_QH,%s_ratio_EK_vs_QH," %(taxa,taxa,taxa,taxa,taxa))
fileOUT_EK_QH.write("\n")


## 3.5 ## PROT FYMINK_GARP
fileOUT_FYMINK_GARP=open("17_FYMINK_GARP.csv","w")
fileOUT_FYMINK_GARP.write("LOCUS,")
for taxa in LT:
    fileOUT_FYMINK_GARP.write("%s_count_FYMINK,%s_prop_FYMINK,%s_count_GARP,%s_prop_GARP," %(taxa,taxa,taxa,taxa))
fileOUT_FYMINK_GARP.write("\n")


## 3.6 ## PROT AVLIMFYW
fileOUT_AVLIMFYW=open("18_AVLIMFYW.csv","w")
fileOUT_AVLIMFYW.write("LOCUS,")
for taxa in LT:
    fileOUT_AVLIMFYW.write("%s_count_AVLIMFYW,%s_prop_AVLIMFYW,%s_count_AVLIM,%s_prop_AVLIM,%s_count_FYW,%s_prop_FYW," %(taxa,taxa,taxa,taxa,taxa,taxa))
fileOUT_AVLIMFYW.write("\n")

## 3.7 ## PROT STNQ
fileOUT_STNQ=open("19_STNQ.csv","w")
fileOUT_STNQ.write("LOCUS,")
for taxa in LT:
    fileOUT_STNQ.write("%s_count_STNQ,%s_prop_STNQ," %(taxa,taxa))
fileOUT_STNQ.write("\n")

## 3.8 ## PROT RHKDE
fileOUT_RHKDE=open("20_RHKDE.csv","w")
fileOUT_RHKDE.write("LOCUS,")
for taxa in LT:
    fileOUT_RHKDE.write("%s_count_RHKDE,%s_prop_RHKDE,%s_count_RHK,%s_prop_RHK,%s_count_DE,%s_prop_DE," %(taxa,taxa,taxa,taxa,taxa,taxa))
fileOUT_RHKDE.write("\n")

## 3.9 ## PROT DIDER CRITERIA
fileOUT_DIDIER=open("21_DIDIER_CRITERIA.csv","w")
fileOUT_DIDIER.write("LOCUS,")
for taxa in LT:
    fileOUT_DIDIER.write("%s_count_APGC,%s_prop_APGC,%s_count_AC,%s_prop_AC,%s_count_VLIM,%s_prop_VLIM,%s_ratio_APGC_vs_VLIM,%s_ratio_AC_vs_VLIM," %(taxa,taxa,taxa,taxa,taxa,taxa,taxa,taxa))
fileOUT_DIDIER.write("\n")

## 3.10 ## PROT Total residue weight
fileOUT_TotalResidueWeight=open("22_TotalResidueWeight.csv","w")
fileOUT_TotalResidueWeight.write("LOCUS,")
for taxa in LT:
    fileOUT_TotalResidueWeight.write("%s_Total_Residue_Weight," %taxa)
fileOUT_TotalResidueWeight.write("\n")

## 3.11 ## PROT Total residue volume
fileOUT_TotalResidueVolume=open("23_TotalResidueVolume.csv","w")
fileOUT_TotalResidueVolume.write("LOCUS,")
for taxa in LT:
    fileOUT_TotalResidueVolume.write("%s_Total_Residue_Volume," %taxa)
fileOUT_TotalResidueVolume.write("\n")

## 3.12 ## PROT Total partial specific volume
fileOUT_TotalPartialSpecificVolume=open("24_TotalPartialSpecificVolume.csv","w")
fileOUT_TotalPartialSpecificVolume.write("LOCUS,")
for taxa in LT:
    fileOUT_TotalPartialSpecificVolume.write("%s_Total_Partial_Specific_Volume," %taxa)
fileOUT_TotalPartialSpecificVolume.write("\n")

## 3.13 ## PROT Total hydratation
fileOUT_TotalHydratation=open("25_TotalHydratation.csv","w")
fileOUT_TotalHydratation.write("LOCUS,")
for taxa in LT:
    fileOUT_TotalHydratation.write("%s_Total_Hydratation," %taxa)
fileOUT_TotalHydratation.write("\n")


#####################
## 4 ## Process Loci
#####################
bash_aa_properties = aa_properties(fileIN_properties)

for locus in Lloci_AA:
    print locus
    path_locus = "%s/%s" %(Path_IN_loci_AA, locus)
    bash = dico(path_locus)

    #print bash
    
    fileOUT_PROT_ALL.write("%s," %locus)
    fileOUT_IVYWREL.write("%s," %locus)
    fileOUT_ERK_DNQTSHA.write("%s," %locus)
    fileOUT_EK_QH.write("%s," %locus)
    fileOUT_FYMINK_GARP.write("%s," %locus)
    fileOUT_AVLIMFYW.write("%s," %locus)
    fileOUT_STNQ.write("%s," %locus)
    fileOUT_RHKDE.write("%s," %locus)
    fileOUT_DIDIER.write("%s," %locus)
    fileOUT_TotalResidueWeight.write("%s," %locus)
    fileOUT_TotalResidueVolume.write("%s," %locus)
    fileOUT_TotalPartialSpecificVolume.write("%s," %locus)
    fileOUT_TotalHydratation.write("%s," %locus)
    
    for taxa in LT:
        seq = bash[taxa]
        prop_K,prop_R,prop_A,prop_F,prop_I,prop_L,prop_M,prop_V,prop_W,prop_N,prop_Q,prop_S,prop_T,prop_H,prop_Y,prop_C,prop_D,prop_E,prop_P,prop_G = aa_composition1(seq)   ### DEF3 ###
        count_IVYWREL,prop_IVYWREL,count_ERK,prop_ERK,count_DNQTSHA,prop_DNQTSHA,ratio_ERK_vs_DNQTSHA,count_EK,prop_EK,count_QH,prop_QH,ratio_EK_vs_QH,count_FYMINK,prop_FYMINK,count_GARP,prop_GARP,count_AVLIMFYW,prop_AVLIMFYW,count_AVLIM,prop_AVLIM,count_FYW,prop_FYW,count_STNQ,prop_STNQ, count_VLIM,prop_VLIM, count_APGC,prop_APGC, count_AC,prop_AC, ratio_APGC_vs_VLIM, ratio_AC_vs_VLIM,count_RHKDE,prop_RHKDE,count_RHK,prop_RHK,count_DE,prop_DE = aa_composition2(seq)   ### DEF4 ###
        Total_Residue_Weight,Total_Residue_Volume,Total_Partial_Specific_Volume,Total_Hydration = sequence_properties_from_aa_properties(seq, bash_aa_properties)   ### DEF6 ###
        
        fileOUT_PROT_ALL.write("%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f," %(prop_K,prop_R,prop_A,prop_F,prop_I,prop_L,prop_M,prop_V,prop_W,prop_N,prop_Q,prop_S,prop_T,prop_H,prop_Y,prop_C,prop_D,prop_E,prop_P,prop_G))
        fileOUT_IVYWREL.write("%.5f,%.5f," %(count_IVYWREL, prop_IVYWREL))
        fileOUT_ERK_DNQTSHA.write("%.5f,%.5f,%.5f,%.5f,%.5f," %(count_ERK,prop_ERK,count_DNQTSHA,prop_DNQTSHA,ratio_ERK_vs_DNQTSHA))
        fileOUT_EK_QH.write("%.5f,%.5f,%.5f,%.5f,%.5f," %(count_EK,prop_EK,count_QH,prop_QH,ratio_EK_vs_QH))
        fileOUT_FYMINK_GARP.write("%.5f,%.5f,%.5f,%.5f," %(count_FYMINK,prop_FYMINK,count_GARP,prop_GARP))
        fileOUT_AVLIMFYW.write("%.5f,%.5f,%.5f,%.5f,%.5f,%.5f," %(count_AVLIMFYW,prop_AVLIMFYW,count_AVLIM,prop_AVLIM,count_FYW,prop_FYW))
        fileOUT_STNQ.write("%.5f,%.5f," %(count_STNQ,prop_STNQ))
        fileOUT_RHKDE.write("%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,"%(count_RHKDE,prop_RHKDE,count_RHK,prop_RHK,count_DE,prop_DE))
        fileOUT_DIDIER.write("%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f," %(count_APGC,prop_APGC,count_AC,prop_AC,count_VLIM,prop_VLIM,ratio_APGC_vs_VLIM,ratio_AC_vs_VLIM))
        fileOUT_TotalResidueWeight.write("%.5f," %Total_Residue_Weight)
        fileOUT_TotalResidueVolume.write("%.5f," %Total_Residue_Volume)
        fileOUT_TotalPartialSpecificVolume.write("%.5f," %(Total_Partial_Specific_Volume))
        fileOUT_TotalHydratation.write("%.5f," % Total_Hydration)
        
    ## END LINE
    fileOUT_PROT_ALL.write("\n")
    fileOUT_IVYWREL.write("\n")
    fileOUT_ERK_DNQTSHA.write("\n")
    fileOUT_EK_QH.write("\n")
    fileOUT_FYMINK_GARP.write("\n")
    fileOUT_AVLIMFYW.write("\n")
    fileOUT_STNQ.write("\n")
    fileOUT_RHKDE.write("\n")    
    fileOUT_DIDIER.write("\n")
    fileOUT_TotalResidueWeight.write("\n")
    fileOUT_TotalResidueVolume.write("\n")
    fileOUT_TotalPartialSpecificVolume.write("\n")
    fileOUT_TotalHydratation.write("\n")







#     # Proportions AA
#     ApS_prop_K,ApS_prop_R,ApS_prop_A,ApS_prop_F,ApS_prop_I,ApS_prop_L,ApS_prop_M,ApS_prop_V,ApS_prop_W,ApS_prop_N,ApS_prop_Q,ApS_prop_S,ApS_prop_T,ApS_prop_H,ApS_prop_Y,ApS_prop_C,ApS_prop_D,ApS_prop_E,ApS_prop_P,ApS_prop_G = aa_composition1(ApS_seq)   ### DEF3 ###
#     # Criterias for testing thermal adaptation
#     ApS_count_IVYWREL,ApS_prop_IVYWREL,ApS_count_ERK,ApS_prop_ERK,ApS_count_DNQTSHA,ApS_prop_DNQTSHA,ApS_ratio_ERK_vs_DNQTSHA,ApS_count_EK,ApS_prop_EK,ApS_count_QH,ApS_prop_QH,ApS_ratio_EK_vs_QH,ApS_count_FYMINK,ApS_prop_FYMINK,ApS_count_GARP,ApS_prop_GARP,ApS_prop_AVLIFYW,ApS_prop_AVLI,ApS_prop_FYW,ApS_prop_STNQ = aa_composition2(ApS_seq)   ### DEF4 ###
    
#     ## 3.2. Get AA composition for Ps
#     # Proportions AA
#     Ps_prop_K,Ps_prop_R,Ps_prop_A,Ps_prop_F,Ps_prop_I,Ps_prop_L,Ps_prop_M,Ps_prop_V,Ps_prop_W,Ps_prop_N,Ps_prop_Q,Ps_prop_S,Ps_prop_T,Ps_prop_H,Ps_prop_Y,Ps_prop_C,Ps_prop_D,Ps_prop_E,Ps_prop_P,Ps_prop_G = aa_composition1(Ps_seq)   ### DEF3 ###
#     # Criterias for testing thermal adaptation
#     Ps_count_IVYWREL,Ps_prop_IVYWREL,Ps_count_ERK,Ps_prop_ERK,Ps_count_DNQTSHA,Ps_prop_DNQTSHA,Ps_ratio_ERK_vs_DNQTSHA,Ps_count_EK,Ps_prop_EK,Ps_count_QH,Ps_prop_QH,Ps_ratio_EK_vs_QH,Ps_count_FYMINK,Ps_prop_FYMINK,Ps_count_GARP,Ps_prop_GARP,Ps_prop_AVLIFYW,Ps_prop_AVLI,Ps_prop_FYW,Ps_prop_STNQ = aa_composition2(Ps_seq)   ### DEF4 ###

#     ## 3.3. Get AA composition for Pg
#     # Proportions AA
#     Pg_prop_K,Pg_prop_R,Pg_prop_A,Pg_prop_F,Pg_prop_I,Pg_prop_L,Pg_prop_M,Pg_prop_V,Pg_prop_W,Pg_prop_N,Pg_prop_Q,Pg_prop_S,Pg_prop_T,Pg_prop_H,Pg_prop_Y,Pg_prop_C,Pg_prop_D,Pg_prop_E,Pg_prop_P,Pg_prop_G = aa_composition1(Pg_seq)   ### DEF3 ###
#     # Criterias for testing thermal adaptation
#     Pg_count_IVYWREL,Pg_prop_IVYWREL,Pg_count_ERK,Pg_prop_ERK,Pg_count_DNQTSHA,Pg_prop_DNQTSHA,Pg_ratio_ERK_vs_DNQTSHA,Pg_count_EK,Pg_prop_EK,Pg_count_QH,Pg_prop_QH,Pg_ratio_EK_vs_QH,Pg_count_FYMINK,Pg_prop_FYMINK,Pg_count_GARP,Pg_prop_GARP,Pg_prop_AVLIFYW,Pg_prop_AVLI,Pg_prop_FYW,Pg_prop_STNQ=aa_composition2(Pg_seq)  ### DEF4 ### 


#     ## 3.2. Print OUTPUT
#     # Proportions AA
#     fileOUT_PROT_ALL.write("%s,%s,%s,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n"  %(ApS_name, Ps_name, Pg_name, ApS_prop_K,ApS_prop_R,ApS_prop_A,ApS_prop_F,ApS_prop_I,ApS_prop_L,ApS_prop_M,ApS_prop_V,ApS_prop_W,ApS_prop_N,ApS_prop_Q,ApS_prop_S,ApS_prop_T,ApS_prop_H,ApS_prop_Y,ApS_prop_C,ApS_prop_D,ApS_prop_E,ApS_prop_P,ApS_prop_G,Ps_prop_K,Ps_prop_R,Ps_prop_A,Ps_prop_F,Ps_prop_I,Ps_prop_L,Ps_prop_M,Ps_prop_V,Ps_prop_W,Ps_prop_N,Ps_prop_Q,Ps_prop_S,Ps_prop_T,Ps_prop_H,Ps_prop_Y,Ps_prop_C,Ps_prop_D,Ps_prop_E,Ps_prop_P,Ps_prop_G,Pg_prop_K,Pg_prop_R,Pg_prop_A,Pg_prop_F,Pg_prop_I,Pg_prop_L,Pg_prop_M,Pg_prop_V,Pg_prop_W,Pg_prop_N,Pg_prop_Q,Pg_prop_S,Pg_prop_T,Pg_prop_H,Pg_prop_Y,Pg_prop_C,Pg_prop_D,Pg_prop_E,Pg_prop_P,Pg_prop_G))

#     # Criterias for testing thermal adaptation    
#     fileOUT_PROT_Criteria.write("%s,%s,%s,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n" %(ApS_name,Ps_name,Pg_name,ApS_count_IVYWREL,ApS_prop_IVYWREL,ApS_count_ERK,ApS_prop_ERK,ApS_count_DNQTSHA,ApS_prop_DNQTSHA,ApS_ratio_ERK_vs_DNQTSHA,ApS_count_EK,ApS_prop_EK,ApS_count_QH,ApS_prop_QH,ApS_ratio_EK_vs_QH,ApS_count_FYMINK,ApS_prop_FYMINK,ApS_count_GARP,ApS_prop_GARP,Ps_count_IVYWREL,Ps_prop_IVYWREL,Ps_count_ERK,Ps_prop_ERK,Ps_count_DNQTSHA,Ps_prop_DNQTSHA,Ps_ratio_ERK_vs_DNQTSHA,Ps_count_EK,Ps_prop_EK,Ps_count_QH,Ps_prop_QH,Ps_ratio_EK_vs_QH,Ps_count_FYMINK,Ps_prop_FYMINK,Ps_count_GARP,Ps_prop_GARP,Pg_count_IVYWREL,Pg_prop_IVYWREL,Pg_count_ERK,Pg_prop_ERK,Pg_count_DNQTSHA,Pg_prop_DNQTSHA,Pg_ratio_ERK_vs_DNQTSHA,Pg_count_EK,Pg_prop_EK,Pg_count_QH,Pg_prop_QH,Pg_ratio_EK_vs_QH,Pg_count_FYMINK,Pg_prop_FYMINK,Pg_count_GARP,Pg_prop_GARP,ApS_prop_AVLIFYW,ApS_prop_AVLI,ApS_prop_FYW,ApS_prop_STNQ,Ps_prop_AVLIFYW,Ps_prop_AVLI,Ps_prop_FYW,Ps_prop_STNQ,Pg_prop_AVLIFYW,Pg_prop_AVLI,Pg_prop_FYW,Pg_prop_STNQ))


# ## 4 ## Get sequence properties based on aa properties

# bash_aa_properties = aa_properties(fileIN_properties)                                         ### DEF5 ###


# for sublist in trio_list:
    
#     ApS_name = sublist[0]
#     ApS_seq = sublist[1]
#     Ps_name = sublist[2]
#     Ps_seq = sublist[3]
#     Pg_name = sublist[4]
#     Pg_seq = sublist[5]

    
#     ApS_Total_Residue_Weight,ApS_Total_Residue_Volume,ApS_Total_Partial_specific_volume,ApS_Total_Hydration =  sequence_properties_from_aa_properties(ApS_seq, bash_aa_properties)                  ### DEF6 ###
#     Ps_Total_Residue_Weight,Ps_Total_Residue_Volume,Ps_Total_Partial_specific_volume,Ps_Total_Hydration =  sequence_properties_from_aa_properties(Ps_seq, bash_aa_properties)                       ### DEF6 ###
#     Pg_Total_Residue_Weight,Pg_Total_Residue_Volume,Pg_Total_Partial_specific_volume,Pg_Total_Hydration =  sequence_properties_from_aa_properties(Pg_seq, bash_aa_properties)                       ### DEF6 ###

#     fileOUT_PROT_Properties.write("%s,%s,%s,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n" %(ApS_name, Ps_name, Pg_name, ApS_Total_Residue_Weight, Ps_Total_Residue_Weight, Pg_Total_Residue_Weight, ApS_Total_Residue_Volume, Ps_Total_Residue_Volume, Pg_Total_Residue_Volume, ApS_Total_Partial_specific_volume, Ps_Total_Partial_specific_volume, Pg_Total_Partial_specific_volume, ApS_Total_Hydration, Ps_Total_Hydration, Pg_Total_Hydration))

    
    

