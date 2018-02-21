#!/usr/bin/env python
# coding:utf-8

# Command line : ./S02a_codon_counting.py <concatenation_file_from_phylogeny>

#python codoncounting.py alignement_file
#Warning!! alignement_file should always be preceded by a path, at least ./
#Warning!! pairs.txt should be located in the same folder as alignement_file

#the script counts the number of codons, amino acids, and types of amino acids in sequences, as well as the mutation bias from one item to another between 2 sequences
#counting is then compared to empirical p-values, obtained from bootstrapped sequences obtained from a subset of sequences
     #the script takes as input the DNA alignment (fasta format): python codoncounting.py file_path.fasta
     #in the output files, the pvalues indicate the position of the observed data in a distribution of empirical countings obtained from a resample of the data. Values above 0.95 indicate a significantly higher counting, values under 0.05 a significantly lower counting

#the script automatically reads the sequences to compare from a file that must be called pairs.txt and located with the .fasta file
#in the pairs.txt file, sequences (let's assume X, Y, Z, U, V) pairs must be written as 'X Y\nU V\nZ V'
#in this case, codoncounting will count the occurence of codons, amino acids, and types of amino acids in X, U, Z, and count the mutation bias from Y to X, V to U and V to Z
#you can add comments in the pairs.txt file inbetween lines, beginning with '#'. E.G. 'X Y\n#This is my comment\nU V\nZ V'
#X, Y, Z, U, V must be character strings contained in the sequences names in the .fasta file (and be specific to each of them)
     #in pairs.txt, you must write how should be built the bootstrapped resampling of sequences. This must be formated as:'X Y\nbackground: length iterration plusminus listofspecies\nU V\nZ V', explanation below
     #backgrounds must be excplicitely written in the pairs.txt file (the script still integers default parameters). This implies that the first line of pairs.txt should be a background line
     #by default, once the background has been determined, it will be applied to each subsequent analysis until another background is written
     #e.g. 'background: length1 iterration1 plusminus1 listofspecies1\nU V\nZ V\nbackground: length2 iterration2 plusminus2 listofspecies2\nX Y' the first background is applied to U V and Z V and the 2nd background to X Y


#the script resamples random pairs of aligned codon to determine what countings can be expected under the hypothesis of an homogenous dataset
#countings are performed on each generated random alignement, thousands of alignments allow to draw a gaussian distribution of the countings
#then the script simply checks whether the observed data are within the 5% lowest or 5% highest values of the distribution
     #in background: length iterration plusminus listofspecies
     #-> length is the number of pairs of codons in each generated alignments (effect on the robustness on the countings performed on this alignement)
     #-> iterration is the number of alignments that will be generated (effect on the resolution of the gaussian distribution)
     #-> plusminus can be either '+' or '-', '+' indicates that the following species only must be resampled, '-' that the following species must be excluded from the resampling
     #-> listofspecies is the list of species (names contained in the sequences names from the fasta file) that must be included or excluded from the sampling. You can also write 'all' to include every species (in this case, plusminus parameter is ignored)
     #full example: background 5000 10000 + melanogaster elegans sapiens
#iterration shouldn't be lower that 1000 to have a relatively smooth gaussian distribution, length shouldn't be lower as 1000 to detect codons with relatively low occurence (<1%)
#for the list of species, you can try to form subgroups depending on the studied parameter (e.g. comparing a terrestrial species with a background composed of marine species)


#added: also counts GC3 and IVYWREL
#added: also counts GC12, CvP bias and EK/QH
#added: also counts purine load and PAYRE/SDGM


def viable(seqs,pos):
  r=0
  compt=0
  for i in range(len(seqs)):
    if i%2==1:
      if not '-' in seqs[i][pos:pos+3]:
        compt+=1
    if compt>1:
      r=1
  return r

def substrcountings(table,pvalues):

  string=[]
  names=[]
  numbers=[]
  stats=['pvalues']
  for f in table:
    names.append(',%s' % f)
    numbers.append(',%f' % table[f])
    stats.append(',%f' % pvalues[f])
  names.append('\n')
  numbers.append('\n')
  stats.append('\n')
  string=names+numbers+stats  

  return string

def substrbiases(table,pvalues):

  string=[]
  names=[]
  numbers=[]
  stats=['pvalues']
  for f in table:
    for g in table[f]:
      names.append(',%s>%s' % (g,f))
      numbers.append(',%f' % table[f][g])
      stats.append(',%f' % pvalues[f][g])
  names.append('\n')
  numbers.append('\n')
  stats.append('\n')  
  string=names+numbers+stats  

  return string

def strcountings(codons, aa, classif,codonspvalues,aapvalues,classifpvalues,GC3,GC12,IVYWREL,EKQH,PAYRESDGM,purineload,CvP):

  subcodons=['occurence of codons\n']+substrcountings(codons,codonspvalues)
  subaa=subcodons+['\noccurence of amino_acids\n']+substrcountings(aa,aapvalues)
  subclassif=subaa+['\noccurence of amino_acid types\n']+substrcountings(classif,classifpvalues)

  string=subclassif+[('\nGC3	GC12	IVYWREL	EK/QH	PAYRE/SDGM	Purine_load	CvP\n%f,%f,%f,%f,%f,%f,%f\n' % (GC3, GC12, IVYWREL, EKQH, PAYRESDGM, purineload, CvP))]

  return string


def strbiases(codons, aa, classif,codonspvalues,aapvalues,classifpvalues,fileaa, fileaatype, p):

  # Writing to files (tab separated; the tab characters are in the lists values)

  aafreq, aapval = substrbiases(aa,aapvalues)[401:801], substrbiases(aa,aapvalues)[803:-1]
  aatypefreq, aatypepval = substrbiases(classif,classifpvalues)[17:33], substrbiases(classif,classifpvalues)[35:-1]
      
  fileaa.write("%s>%s" %(p[1], p[0]))
  for value in aafreq:
    fileaa.write(str(value))
  fileaa.write("\n")
  fileaa.write("%s>%s_pvalue" %(p[1], p[0]))
  for value in aapval:
    fileaa.write(str(value))
  fileaa.write("\n")
  
  fileaatype.write("%s>%s" %(p[1], p[0]))
  for value in aatypefreq:
    fileaatype.write(str(value))
  fileaatype.write("\n")
  fileaatype.write("%s>%s_pvalue" %(p[1], p[0]))
  for value in aatypepval:
    fileaatype.write(str(value))
  fileaatype.write("\n")

  # Origin code

  subcodons=['codons mutations biases\n']+substrbiases(codons,codonspvalues)  
  subaa=subcodons+['\namino acids mutation biases\n']+substrbiases(aa,aapvalues)
  subclassif=subaa+['\ntypes of amino_acids mutation biasecodoncountingv22.pys\n']+substrbiases(classif,classifpvalues)  

  return subclassif

def testpvalue(bootstrap,value,iterration): #computes where the observed value is located in the expected counting distribution

  maxval=iterration-1
  minval=0
  testval=(maxval+minval)/2
  print bootstrap[testval]
  while maxval-minval > 1:
    if value > bootstrap[testval]:
      minval=testval
      testval=(maxval+minval)/2
    elif value < bootstrap[testval]:
      maxval=testval
      testval=(maxval+minval)/2
    else:
      break
  pvalue=(testval+1)/float(iterration+1)

  return pvalue

def gettables(content,reversecode,code,classif): #generates the tables contening all the countings

  if content==[]:
    codonscount={k:[] for k in reversecode}
    aacount={k:[] for k in code}
    aaclassifcount={k:[] for k in classif}
  
    codons={}
    for k in reversecode:
      codons[k]={}
      for q in reversecode:
        codons[k][q]=[]
    aa={}
    for k in code:
      aa[k]={}
      for q in code:
        aa[k][q]=[]
    aaclassif={}
    for k in classif:
      aaclassif[k]={}
      for q in classif:
        aaclassif[k][q]=[]

  elif content==0:
    codonscount={k:0 for k in reversecode}
    aacount={k:0 for k in code}
    aaclassifcount={k:0 for k in classif}
  
    codons={}
    for k in reversecode:
      codons[k]={}
      for q in reversecode:
        codons[k][q]=0
    aa={}
    for k in code:
      aa[k]={}
      for q in code:
        aa[k][q]=0
    aaclassif={}
    for k in classif:
      aaclassif[k]={}
      for q in classif:
        aaclassif[k][q]=0

  return codonscount, aacount, aaclassifcount, codons, aa, aaclassif

def countings(seq1,seq2,code,classif,reversecode,reverseclassif):
# COMPTAGES (output counts.txt) functions used : gettables()
#countings actually counts occurence and mutation bias of codons, amino acids and types of amino acids
  
  codonscount, aacount, aaclassifcount, codons, aa, aaclassif=gettables(0,reversecode,code,classif)  
  codonscount2, aacount2, aaclassifcount2, _, _, _=gettables(0,reversecode,code,classif)
 
  G12=0
  C12=0
  G3=0
  C3=0
  A=0
  T=0

  i=0
  compcodons=0
  while i<len(seq1)-1:
    if not '-' in seq1[i:i+3]: #coutings of occurences is obtained from the maximum of full codons being in the sequence
      compcodons+=1

      if (seq1[i]=='g'):
        G12+=1
      if (seq1[i]=='c'):
        C12+=1
      if (seq1[i+1]=='g'):
        G12+=1
      if (seq1[i+1]=='c'):
        C12+=1
      if (seq1[i+2]=='g'):
        G3+=1
      if (seq1[i+2]=='c'):
        C3+=1

      if (seq1[i]=='a'):
        A+=1
      if (seq1[i]=='t'):
        T+=1
      if (seq1[i+1]=='a'):
        A+=1
      if (seq1[i+1]=='t'):
        T+=1
      if (seq1[i+2]=='a'):
        A+=1
      if (seq1[i+2]=='t'):
        T+=1

      codonscount[seq1[i:i+3]]+=1
      aacount[reversecode[seq1[i:i+3]]]+=1
      aaclassifcount[reverseclassif[reversecode[seq1[i:i+3]]]]+=1
    if (not '-' in seq1[i:i+3]) and (not '-' in seq2[i:i+3]) and (not seq1[i:i+3]==seq2[i:i+3]): #mutation biases are count from non ambiguous pairs of codons in the 2 sequences
      codons[seq1[i:i+3]][seq2[i:i+3]]+=1
      codonscount2[seq2[i:i+3]]+=1
      aacount2[reversecode[seq2[i:i+3]]]+=1
      aaclassifcount2[reverseclassif[reversecode[seq2[i:i+3]]]]+=1
    i+=3
  
  IVYWREL=aacount['ile']+aacount['val']+aacount['tyr']+aacount['trp']+aacount['arg']+aacount['glu']+aacount['leu']
  EKQH=(aacount['glu']+aacount['lys'])/max([float(aacount['gln']+aacount['his']),1])
  PAYRESDGM=(aacount['pro']+aacount['ala']+aacount['tyr']+aacount['arg']+aacount['glu'])/max([float(aacount['ser']+aacount['asp']+aacount['gly']+aacount['met']),1])
  
  for i in codons:
    for j in codons[i]:
      if not reversecode[i]==reversecode[j]:
        aa[reversecode[i]][reversecode[j]]+=codons[i][j]

  for i in aa:
    for j in aa[i]:
      if not reverseclassif[i]==reverseclassif[j]:
        aaclassif[reverseclassif[i]][reverseclassif[j]]+=aa[i][j]
        

  for f in codons: #mutations from codon C in sequence X to codon c in sequence Y are normalized by the total number of codons C in sequence X
    for g in codons[f]:
      codons[f][g]=codons[f][g]/float(max([codonscount2[g],1]))
  for f in aa:
    for g in aa[f]:
      aa[f][g]=aa[f][g]/float(max([aacount2[g],1]))
  for f in aaclassif:
    for g in aaclassif[f]:
      aaclassif[f][g]=aaclassif[f][g]/float(max([aaclassifcount2[g],1]))
  for f in codonscount: #occurences are normalized by the total number of non ambiguous codons in the sequence
    codonscount[f]=codonscount[f]/float(compcodons)
  for f in aacount:
    aacount[f]=aacount[f]/float(compcodons)
  for f in aaclassifcount:
    aaclassifcount[f]=aaclassifcount[f]/float(compcodons)
    
  for f in codons:
    for g in codons[f]:
        if f==g:
          break
        if (codons[g][f]==0) and (codons[f][g])>0:
          codons[f][g]=1
        elif (codons[g][f]==0) and (codons[f][g])==0:
          codons[f][g]=0
        else:
          x=codons[f][g]/codons[g][f]
          codons[f][g]=-pow(2,1-x)+1

  for f in aa:
    for g in aa[f]:
        if f==g:
          break
        if (aa[g][f]==0) and (aa[f][g])>0:
          aa[f][g]=1
        elif (aa[g][f]==0) and (aa[f][g])==0:
          aa[f][g]=0      
        else:
          x=aa[f][g]/aa[g][f]
          aa[f][g]=-pow(2,1-x)+1

  for f in aaclassif:
    for g in aaclassif[f]:
        if f==g:
          break
        if (aaclassif[g][f]==0) and (aaclassif[f][g])>0:
          aaclassif[f][g]=1
        elif (aaclassif[g][f]==0) and (aaclassif[f][g])==0:
          aaclassif[f][g]=0
        else:
          x=aaclassif[f][g]/aaclassif[g][f]
          aaclassif[f][g]=-pow(2,1-x)+1

  GC3=(G3+C3)/float(compcodons)
  GC12=(G12+C12)/float(2*compcodons)
  purineload=(1000*(G12+G3+A-C12-C3-T))/float(3*compcodons)
  IVYWREL=IVYWREL/float(compcodons)
  CvP=aaclassifcount['charged']-aaclassifcount['polar']
  
  return codonscount, aacount, aaclassifcount, codons, aa, aaclassif, GC3, GC12, IVYWREL, EKQH, PAYRESDGM, purineload, CvP

def sampling(inputfile,length,iterration,plusminus,species,code,classif,reversecode,reverseclassif):
# functions used : gettables, countings
#sampling provides 'iterations' pairs of sequences of 'length' non ambiguous codons obtained from the dataset specified by 'species' and 'plusminus'
#sort of bootstrap

  codonscount, aacount, aaclassifcount, codons, aa, aaclassif=gettables([],reversecode,code,classif)

  inputfile.seek(0) #generates the bootstrapped sequences
  lines=[]
  comp=-1
  if species==['all']:
    while 1:
      line=inputfile.readline()[:-1]
      if not line:
        break
      lines.append(line)
      comp+=1
  else:
    if plusminus=='+':
       while 1:
        line1=inputfile.readline()[:-1]
        if not line1:
          break
        line2=inputfile.readline()[:-1]
        if any(spe in line1 for spe in species):
          lines.append(line1)
          lines.append(line2)
          comp+=2
    elif plusminus=='-':
       while 1:
        line1=inputfile.readline()[:-1]
        if not line1:
          break
        line2=inputfile.readline()[:-1]
        if not any(spe in line1 for spe in species):
          lines.append(line1)
          lines.append(line2)
          comp+=2
  l=len(lines[1])-1
  
  for z in range(iterration):
    if z%(iterration/10)==0:
      print str(10*z/(iterration/10))+' %'
    seqa=[]
    seqb=[]
    for i in range(length):
      site=random.randrange(0,l,3)
      while viable(lines,site)==0:
        site=random.randrange(0,l,3)
      a='---'
      b='---'
      while ('-' in a) or ('-' in b):
        posa=2
        posb=2
        while posa==posb:
          posa=random.randrange(1,comp+1,2)
          posb=random.randrange(1,comp+1,2)
        a=lines[posa][site:site+3]
        b=lines[posb][site:site+3]
      seqa.append(a)
      seqb.append(b)
    seqa=(''.join(seqa)).lower()
    seqb=(''.join(seqb)).lower()

    codonscounttemp, aacounttemp, aaclassifcounttemp, codonstemp, aatemp, aaclassiftemp, _, _, _, _, _, _, _=countings(seqa,seqb,code,classif,reversecode,reverseclassif)
    for f in codonscount:
      codonscount[f].append(codonscounttemp[f])
    for f in aacount:
      aacount[f].append(aacounttemp[f])
    for f in aaclassifcount:
      aaclassifcount[f].append(aaclassifcounttemp[f])
    for f in codons:
      for g in codons[f]:
        codons[f][g].append(codonstemp[f][g])
    for f in aa:
      for g in aa[f]:
        aa[f][g].append(aatemp[f][g])
    for f in aaclassif:
      for g in aaclassif[f]:
        aaclassif[f][g].append(aaclassiftemp[f][g]) #counts the occurences and mutation biases for each the bootstrapped sequence

  print '100 %'

  for f in codonscount:
    codonscount[f].sort()
  for f in aacount:
    aacount[f].sort()
  for f in aaclassifcount:
    aaclassifcount[f].sort()
  for f in codons:
    for g in codons[f]:
      codons[f][g].sort()
  for f in aa:
    for g in aa[f]:
      aa[f][g].sort()
  for f in aaclassif:
    for g in aaclassif[f]:
      aaclassif[f][g].sort()
    
  return codonscount, aacount, aaclassifcount, codons, aa, aaclassif
 
################
######RUN#######
################
 
import string, os, sys, re, random
from math import pow

codons_counts = open("codons_counts.csv", "w")
aa_counts = open("aa_counts.csv", "w")
aatypes_counts = open("aatypes_counts.csv", "w")
gc_counts = open("gc_counts.csv", "w")

# ou bien utiliser le reversecode plus bas
column_index_codons = "Species,tat,tgt,tct,ttt,tgc,tgg,tac,ttc,tcg,tta,ttg,tcc,tca,gca,gta,gcc,gtc,gcg,gtg,caa,gtt,gct,acc,ggt,cga,cgc,gat,aag,cgg,act,ggg,gga,ggc,gag,aaa,gac,cgt,gaa,ctt,atg,aca,acg,atc,aac,ata,agg,cct,agc,aga,cat,aat,att,ctg,cta,ctc,cac,ccg,agt,cag,cca,ccc\n"
column_index_aa = "Species,cys,asn,his,ile,ser,gln,lys,met,pro,thr,phe,ala,gly,val,leu,asp,arg,trp,glu,tyr\n"
column_index_aatypes = "Species,aromatics,polar,unpolar,charged\n"
column_index_gc = "Species,GC3,GC12,IVYWREL,EKQH,PAYRESDGM,purineload,CvP\n"

codons_counts.write(column_index_codons)
aa_counts.write(column_index_aa)
aatypes_counts.write(column_index_aatypes)
gc_counts.write(column_index_gc)

aa_transitions = open("aa_transitions.csv", "w")
aatypes_transitions = open("aatypes_transitions.csv", "w")

column_index_aa_transitions = "Species,cys>cys,asn>cys,his>cys,ile>cys,ser>cys,gln>cys,lys>cys,met>cys,pro>cys,thr>cys,phe>cys,ala>cys,gly>cys,val>cys,leu>cys,asp>cys,arg>cys,trp>cys,glu>cys,tyr>cys,cys>asn,asn>asn,his>asn,ile>asn,ser>asn,gln>asn,lys>asn,met>asn,pro>asn,thr>asn,phe>asn,ala>asn,gly>asn,val>asn,leu>asn,asp>asn,arg>asn,trp>asn,glu>asn,tyr>asn,cys>his,asn>his,his>his,ile>his,ser>his,gln>his,lys>his,met>his,pro>his,thr>his,phe>his,ala>his,gly>his,val>his,leu>his,asp>his,arg>his,trp>his,glu>his,tyr>his,cys>ile,asn>ile,his>ile,ile>ile,ser>ile,gln>ile,lys>ile,met>ile,pro>ile,thr>ile,phe>ile,ala>ile,gly>ile,val>ile,leu>ile,asp>ile,arg>ile,trp>ile,glu>ile,tyr>ile,cys>ser,asn>ser,his>ser,ile>ser,ser>ser,gln>ser,lys>ser,met>ser,pro>ser,thr>ser,phe>ser,ala>ser,gly>ser,val>ser,leu>ser,asp>ser,arg>ser,trp>ser,glu>ser,tyr>ser,cys>gln,asn>gln,his>gln,ile>gln,ser>gln,gln>gln,lys>gln,met>gln,pro>gln,thr>gln,phe>gln,ala>gln,gly>gln,val>gln,leu>gln,asp>gln,arg>gln,trp>gln,glu>gln,tyr>gln,cys>lys,asn>lys,his>lys,ile>lys,ser>lys,gln>lys,lys>lys,met>lys,pro>lys,thr>lys,phe>lys,ala>lys,gly>lys,val>lys,leu>lys,asp>lys,arg>lys,trp>lys,glu>lys,tyr>lys,cys>met,asn>met,his>met,ile>met,ser>met,gln>met,lys>met,met>met,pro>met,thr>met,phe>met,ala>met,gly>met,val>met,leu>met,asp>met,arg>met,trp>met,glu>met,tyr>met,cys>pro,asn>pro,his>pro,ile>pro,ser>pro,gln>pro,lys>pro,met>pro,pro>pro,thr>pro,phe>pro,ala>pro,gly>pro,val>pro,leu>pro,asp>pro,arg>pro,trp>pro,glu>pro,tyr>pro,cys>thr,asn>thr,his>thr,ile>thr,ser>thr,gln>thr,lys>thr,met>thr,pro>thr,thr>thr,phe>thr,ala>thr,gly>thr,val>thr,leu>thr,asp>thr,arg>thr,trp>thr,glu>thr,tyr>thr,cys>phe,asn>phe,his>phe,ile>phe,ser>phe,gln>phe,lys>phe,met>phe,pro>phe,thr>phe,phe>phe,ala>phe,gly>phe,val>phe,leu>phe,asp>phe,arg>phe,trp>phe,glu>phe,tyr>phe,cys>ala,asn>ala,his>ala,ile>ala,ser>ala,gln>ala,lys>ala,met>ala,pro>ala,thr>ala,phe>ala,ala>ala,gly>ala,val>ala,leu>ala,asp>ala,arg>ala,trp>ala,glu>ala,tyr>ala,cys>gly,asn>gly,his>gly,ile>gly,ser>gly,gln>gly,lys>gly,met>gly,pro>gly,thr>gly,phe>gly,ala>gly,gly>gly,val>gly,leu>gly,asp>gly,arg>gly,trp>gly,glu>gly,tyr>gly,cys>val,asn>val,his>val,ile>val,ser>val,gln>val,lys>val,met>val,pro>val,thr>val,phe>val,ala>val,gly>val,val>val,leu>val,asp>val,arg>val,trp>val,glu>val,tyr>val,cys>leu,asn>leu,his>leu,ile>leu,ser>leu,gln>leu,lys>leu,met>leu,pro>leu,thr>leu,phe>leu,ala>leu,gly>leu,val>leu,leu>leu,asp>leu,arg>leu,trp>leu,glu>leu,tyr>leu,cys>asp,asn>asp,his>asp,ile>asp,ser>asp,gln>asp,lys>asp,met>asp,pro>asp,thr>asp,phe>asp,ala>asp,gly>asp,val>asp,leu>asp,asp>asp,arg>asp,trp>asp,glu>asp,tyr>asp,cys>arg,asn>arg,his>arg,ile>arg,ser>arg,gln>arg,lys>arg,met>arg,pro>arg,thr>arg,phe>arg,ala>arg,gly>arg,val>arg,leu>arg,asp>arg,arg>arg,trp>arg,glu>arg,tyr>arg,cys>trp,asn>trp,his>trp,ile>trp,ser>trp,gln>trp,lys>trp,met>trp,pro>trp,thr>trp,phe>trp,ala>trp,gly>trp,val>trp,leu>trp,asp>trp,arg>trp,trp>trp,glu>trp,tyr>trp,cys>glu,asn>glu,his>glu,ile>glu,ser>glu,gln>glu,lys>glu,met>glu,pro>glu,thr>glu,phe>glu,ala>glu,gly>glu,val>glu,leu>glu,asp>glu,arg>glu,trp>glu,glu>glu,tyr>glu,cys>tyr,asn>tyr,his>tyr,ile>tyr,ser>tyr,gln>tyr,lys>tyr,met>tyr,pro>tyr,thr>tyr,phe>tyr,ala>tyr,gly>tyr,val>tyr,leu>tyr,asp>tyr,arg>tyr,trp>tyr,glu>tyr,tyr>tyr\n"
column_index_aatypes_transitions = "Species,aromatics>aromatics,polar>aromatics,unpolar>aromatics,charged>aromatics,aromatics>polar,polar>polar,unpolar>polar,charged>polar,aromatics>unpolar,polar>unpolar,unpolar>unpolar,charged>unpolar,aromatics>charged,polar>charged,unpolar>charged,charged>charged\n"

aa_transitions.write(column_index_aa_transitions)
aatypes_transitions.write(column_index_aatypes_transitions)

PATH = sys.argv[1]
length=1000
iterration=100
background=0
speciesboot=['all']
stringcounts=[' ']
stringbiases=[' ']
towrite=1

code={'phe':['ttt','ttc'],'leu':['tta','ttg','ctt','ctc','cta','ctg'],'ile':['att','atc','ata'],'met':['atg'],'val':['gtt','gtc','gta','gtg'],'ser':['tct','tcc','tca','tcg','agt','agc'],'pro':['cct','cca','ccg','ccc'],'thr':['act','acc','aca','acg'],'ala':['gct','gcc','gca','gcg'],'tyr':['tat','tac'],'his':['cat','cac'],'gln':['caa','cag'],'asn':['aat','aac'],'lys':['aaa','aag'],'asp':['gat','gac'],'glu':['gaa','gag'],'cys':['tgt','tgc'],'trp':['tgg'],'arg':['cgt','cgc','cga','cgg','aga','agg'],'gly':['ggt','ggc','gga','ggg']}
classif={'unpolar':['gly','ala','val','leu','met','ile'],'polar':['ser','thr','cys','pro','asn','gln'],'charged':['lys','arg','his','asp','glu'],'aromatics':['phe','tyr','trp']}

reversecode={v:k for k in code for v in code[k]}
reverseclassif={v:k for k in classif for v in classif[k]} 

# --------------------------------------------------------------

pairs=open("pairs.txt","r") #finds the pairs to analyze
pairlist=[]
while 1:
  line=pairs.readline()
  if not line:
    #if not lastline.startswith('#'):
      #pairlist[-1][1]=string.split(lastline,' ')[1]
    break
  #lastline=line
  if not line.startswith('#'):
    if not line.startswith('background:'):
      if background==0:
        pairlist.append([string.split(line,' ')[0],string.split(line,' ')[1][:-1]])
      elif background==1:
        pairlist.append([string.split(line,' ')[0],string.split(line,' ')[1][:-1],length,iterration,plusminus,species])
        background=0
    else:
      parameters=string.split(line[:-1],': ')[1]
      listparameters=string.split(parameters,' ')
      length=int(listparameters[0])
      iterration=int(listparameters[1])
      plusminus=listparameters[2]
      species=listparameters[3:]
      background=1

# --------------------------------------------------------------

concat=open(sys.argv[1],"r")

last_iter = 1
for p in pairlist: #pairs analysis  
  concat.seek(0)
  while 1:
    # Everything lowercase
    line=concat.readline()
    seq=concat.readline()
    if p[0] in line:
      seq1=seq[:-1].lower()
    elif p[1] in line:
      seq2=seq[:-1].lower()
    if not line:
      break
  
  if len(p)>2: #bootstrap simulations if the background changed
    # Only len(pairlist[0]) > 2
    length=p[2]
    iterration=p[3]
    plusminus=p[4]
    species=p[5] # ex : 'all'
    applause='background: '+str(length)+' '+str(iterration)+' '+str(plusminus)+' '+' '.join(species)
    print applause
    stringcounts.append(applause+'\n\n\n')
    stringbiases.append(applause+'\n\n\n')
    # variables below are computed on the first iteration and then used for the rest of the loop
    codonscountboot, aacountboot, aaclassifcountboot, codonsboot, aaboot, aaclassifboot=sampling(concat,length,iterration,plusminus,species,code,classif,reversecode,reverseclassif)
  print str(p[0])+' vs '+str(p[1])

  # -----------------------------
  # Mise en memoire des calculs pour les comptages par sps + Ecriture fichiers de sortie
  codonscount, aacount, aaclassifcount, codons, aa, aaclassif, GC3, GC12, IVYWREL, EKQH, PAYRESDGM, purineload, CvP=countings(seq1,seq2,code,classif,reversecode,reverseclassif)
  codonscountpvalue, aacountpvalue, aaclassifcountpvalue, codonspvalue, aapvalue, aaclassifpvalue=gettables(0,reversecode,code,classif)


  for f in codons: #tests the countings and saves the pvalues
    for g in codons[f]:
      codonspvalue[f][g]=testpvalue(codonsboot[f][g],codons[f][g],iterration)
  for f in aa:
    for g in aa[f]:
      aapvalue[f][g]=testpvalue(aaboot[f][g],aa[f][g],iterration)
  for f in aaclassif:
    for g in aaclassif[f]:
      aaclassifpvalue[f][g]=testpvalue(aaclassifboot[f][g],aaclassif[f][g],iterration)


  for substring in stringcounts: #to not write twice the same counting with the same background
    if (("counting of %s" % p[0]) in substring) and (applause in substring):
      towrite=0

  if towrite:
    for f in codonscount:
      codonscountpvalue[f]=testpvalue(codonscountboot[f],codonscount[f],iterration)
    for f in aacount:
      aacountpvalue[f]=testpvalue(aacountboot[f],aacount[f],iterration)
    for f in aaclassifcount:
      aaclassifcountpvalue[f]=testpvalue(aaclassifcountboot[f],aaclassifcount[f],iterration)
   
    ## Writing countings into separated output files ##

    codons_counts.write(p[0] + ",") # Species name
    #codons_counts.write(",") # next cell
    for value in codonscount.values()[0:-1]:
      codons_counts.write(str(value) + ",") # write codons_counts line for the species
    codons_counts.write(str(codonscount.values()[-1]))
    codons_counts.write("\n") # next line

    codons_counts.write(p[0] + "_pvalue,") # Same species, but line for computed pvalues    
    for value in codonscountpvalue.values()[0:-1]:
      codons_counts.write(str(value) + ",") # write pvalues line for the species
    codons_counts.write(str(codonscountpvalue.values()[-1]))
    codons_counts.write("\n")

    aa_counts.write(p[0] + ",") # Same method as above    
    for value in aacount.values()[0:-1]:
      aa_counts.write(str(value) + ",")
    aa_counts.write(str(aacount.values()[-1]))
    aa_counts.write("\n")

    aa_counts.write(p[0] + "_pvalue,")
    for value in aacountpvalue.values()[0:-1]:
      aa_counts.write(str(value) + ",")
    aa_counts.write(str(aacountpvalue.values()[-1]))
    aa_counts.write("\n")

    aatypes_counts.write(p[0] + ",")
    for value in aaclassifcount.values()[0:-1]:
      aatypes_counts.write(str(value) + ",")
    aatypes_counts.write(str(aaclassifcount.values()[-1]))
    aatypes_counts.write("\n")

    aatypes_counts.write(p[0] + "_pvalue,")
    for value in aaclassifcountpvalue.values()[0:-1]:
      aatypes_counts.write(str(value) + ",")
    aatypes_counts.write(str(aaclassifcountpvalue.values()[-1]))
    aatypes_counts.write("\n")

    gc_counts.write(p[0])
    gc_counts.write(",")
    gc_counts.write(str(GC3)+","+str(GC12)+","+str(IVYWREL)+","+str(EKQH)+","+str(PAYRESDGM)+","+str(purineload)+","+str(CvP)+"\n")    

    """ IMPROVMENT BELOW :
    Countings was not done on the last species of the list. It is now compute when the loop reaches the last iteration, thanks to the countings()
    function : it is called with a switch between the arguments 'seq1' and seq2'.    
    """
    
    if last_iter == len(pairlist):
      codonscount2, aacount2, aaclassifcount2, codons2, aa2, aaclassif2, GC3_b, GC12_b, IVYWREL_b, EKQH_b, PAYRESDGM_b, purineload_b, CvP_b=countings(seq2,seq1,code,classif,reversecode,reverseclassif)
      codonscountpvalue2, aacountpvalue2, aaclassifcountpvalue2, codonspvalue2, aapvalue2, aaclassifpvalue2=gettables(0,reversecode,code,classif)

      for f in codons2: #tests the countings and saves the pvalues
        for g in codons2[f]:
          codonspvalue2[f][g]=testpvalue(codonsboot[f][g],codons[f][g],iterration)
      for f in aa2:
        for g in aa2[f]:
          aapvalue2[f][g]=testpvalue(aaboot[f][g],aa[f][g],iterration)
      for f in aaclassif2:
        for g in aaclassif2[f]:
          aaclassifpvalue2[f][g]=testpvalue(aaclassifboot[f][g],aaclassif[f][g],iterration)

      for f in codonscount2:
        codonscountpvalue2[f]=testpvalue(codonscountboot[f],codonscount[f],iterration)
      for f in aacount2:
        aacountpvalue2[f]=testpvalue(aacountboot[f],aacount[f],iterration)
      for f in aaclassifcount2:
        aaclassifcountpvalue2[f]=testpvalue(aaclassifcountboot[f],aaclassifcount[f],iterration)

      # write last line of each countings file

      codons_counts.write(p[1] + ",") # second species of the couple, the last one to be written in the file
      for value in codonscount2.values()[0:-1]:
        codons_counts.write(str(value) + ",")
      codons_counts.write(str(codonscount2.values()[-1]))
      codons_counts.write("\n")

      codons_counts.write(p[1] + "_pvalue,") # Same species, but line for computed pvalues
      for value in codonscountpvalue2.values()[0:-1]:
        codons_counts.write(str(value) + ",") # write pvalues line for the species
      codons_counts.write(str(codonscountpvalue2.values()[-1]))  
      codons_counts.write("\n")

      aa_counts.write(p[1] + ",")
      for value in aacount2.values()[0:-1]:
        aa_counts.write(str(value) + ",")
      aa_counts.write(str(aacount2.values()[-1]))
      aa_counts.write("\n")

      aa_counts.write(p[1] + "_pvalue,")
      for value in aacountpvalue2.values()[0:-1]:
        aa_counts.write(str(value) + ",") # write pvalues line for the species
      aa_counts.write(str(aacountpvalue2.values()[-1]))
      aa_counts.write("\n")

      aatypes_counts.write(p[1] + ",")
      for value in aaclassifcount2.values()[0:-1]:
        aatypes_counts.write(str(value) + ",")
      aatypes_counts.write(str(aaclassifcount2.values()[-1]))
      aatypes_counts.write("\n")

      aatypes_counts.write(p[1] + "_pvalue,")
      for value in aaclassifcountpvalue2.values()[0:-1]:
        aatypes_counts.write(str(value) + ",") # write pvalues line for the species
      aatypes_counts.write(str(aaclassifcountpvalue2.values()[-1]))
      aatypes_counts.write("\n")

      gc_counts.write(p[1])
      gc_counts.write(",")
      gc_counts.write(str(GC3_b)+","+str(GC12_b)+","+str(IVYWREL_b)+","+str(EKQH_b)+","+str(PAYRESDGM_b)+","+str(purineload_b)+","+str(CvP_b))

    # end writing    
    
    stringcounts=stringcounts[:-1]+[''.join([stringcounts[-1],("counting of %s\n\n" % p[0])+''.join(strcountings(codonscount, aacount, aaclassifcount,codonscountpvalue,aacountpvalue,aaclassifcountpvalue,GC3,GC12,IVYWREL,EKQH,PAYRESDGM,purineload,CvP))+'\n\n'])]
  
  else:
    towrite=1
  # -----------------------------

  stringbiases.append("mutation biases from %s to %s\n\n" % (p[1], p[0]))
  stringbiases=stringbiases+strbiases(codons, aa, aaclassif,codonspvalue,aapvalue,aaclassifpvalue,aa_transitions,aatypes_transitions,p)
  stringbiases.append('\n\n') 
  print 'done'
  last_iter +=1

concat.close()

# --------------------------------------------------------------

stringcounts=''.join(stringcounts)
stringbiases=''.join(stringbiases)

# results=open(os.path.dirname(PATH)+"/"+string(os.path.split(PATH)[1],' ')[0]+'_results.txt',"w")
#results=open('./codoncounting_results.txt',"w")
results1=open('./counts.txt', "w")
results1.write("%s" % stringcounts)
results1.close()
results2=open('./biases.txt', "w")
results2.write("%s" % stringbiases)
results2.close()
#results.write("%s" % stringcounts)
#results.write("\n\n")
#results.write("%s" % stringbiases)
#results.close()

codons_counts.close()
aa_counts.close()
aatypes_counts.close()
gc_counts.close()
aa_transitions.close()
aatypes_transitions.close()
