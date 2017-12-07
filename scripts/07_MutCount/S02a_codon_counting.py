#!/usr/bin/python

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
    names.append('	%s' % f)
    numbers.append('	%f' % table[f])
    stats.append('	%f' % pvalues[f])
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
      names.append('	%s>%s' % (g,f))
      numbers.append('	%f' % table[f][g])
      stats.append('	%f' % pvalues[f][g])
  names.append('\n')
  numbers.append('\n')
  stats.append('\n')  
  string=names+numbers+stats  

  return string

def strcountings(codons, aa, classif,codonspvalues,aapvalues,classifpvalues,GC3,GC12,IVYWREL,EKQH,PAYRESDGM,purineload,CvP):

  subcodons=['occurence of codons\n']+substrcountings(codons,codonspvalues)
  subaa=subcodons+['\noccurence of amino_acids\n']+substrcountings(aa,aapvalues)
  subclassif=subaa+['\noccurence of amino_acid types\n']+substrcountings(classif,classifpvalues)

  string=subclassif+[('\nGC3	GC12	IVYWREL	EK/QH	PAYRE/SDGM	Purine_load	CvP\n%f	%f	%f	%f	%f	%f	%f\n' % (GC3, GC12, IVYWREL, EKQH, PAYRESDGM, purineload, CvP))]

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
column_index_codons = "Species\ttat\ttgt\ttct\tttt\ttgc\ttgg\ttac\tttc\ttcg\ttta\tttg\ttcc\ttca\tgca\tgta\tgcc\tgtc\tgcg\tgtg\tcaa\tgtt\tgct\tacc\tggt\tcga\tcgc\tgat\taag\tcgg\tact\tggg\tgga\tggc\tgag\taaa\tgac\tcgt\tgaa\tctt\tatg\taca\tacg\tatc\taac\tata\tagg\tcct\tagc\taga\tcat\taat\tatt\tctg\tcta\tctc\tcac\tccg\tagt\tcag\tcca\tccc\n"
column_index_aa = "Species\tcys\tasn\this\tile\tser\tgln\tlys\tmet\tpro\tthr\tphe\tala\tgly\tval\tleu\tasp\targ\ttrp\tglu\ttyr\n"
column_index_aatypes = "Species\taromatics\tpolar\tunpolar\tcharged\n"
column_index_gc = "Species\tGC3\tGC12\tIVYWREL\tEKQH\tPAYRESDGM\tpurineload\tCvP\n"

codons_counts.write(column_index_codons)
aa_counts.write(column_index_aa)
aatypes_counts.write(column_index_aatypes)
gc_counts.write(column_index_gc)

aa_transitions = open("aa_transitions.csv", "w")
aatypes_transitions = open("aatypes_transitions.csv", "w")

column_index_aa_transitions = "Species\tcys>cys\tasn>cys\this>cys\tile>cys\tser>cys\tgln>cys\tlys>cys\tmet>cys\tpro>cys\tthr>cys\tphe>cys\tala>cys\tgly>cys\tval>cys\tleu>cys\tasp>cys\targ>cys\ttrp>cys\tglu>cys\ttyr>cys\tcys>asn\tasn>asn\this>asn\tile>asn\tser>asn\tgln>asn\tlys>asn\tmet>asn\tpro>asn\tthr>asn\tphe>asn\tala>asn\tgly>asn\tval>asn\tleu>asn\tasp>asn\targ>asn\ttrp>asn\tglu>asn\ttyr>asn\tcys>his\tasn>his\this>his\tile>his\tser>his\tgln>his\tlys>his\tmet>his\tpro>his\tthr>his\tphe>his\tala>his\tgly>his\tval>his\tleu>his\tasp>his\targ>his\ttrp>his\tglu>his\ttyr>his\tcys>ile\tasn>ile\this>ile\tile>ile\tser>ile\tgln>ile\tlys>ile\tmet>ile\tpro>ile\tthr>ile\tphe>ile\tala>ile\tgly>ile\tval>ile\tleu>ile\tasp>ile\targ>ile\ttrp>ile\tglu>ile\ttyr>ile\tcys>ser\tasn>ser\this>ser\tile>ser\tser>ser\tgln>ser\tlys>ser\tmet>ser\tpro>ser\tthr>ser\tphe>ser\tala>ser\tgly>ser\tval>ser\tleu>ser\tasp>ser\targ>ser\ttrp>ser\tglu>ser\ttyr>ser\tcys>gln\tasn>gln\this>gln\tile>gln\tser>gln\tgln>gln\tlys>gln\tmet>gln\tpro>gln\tthr>gln\tphe>gln\tala>gln\tgly>gln\tval>gln\tleu>gln\tasp>gln\targ>gln\ttrp>gln\tglu>gln\ttyr>gln\tcys>lys\tasn>lys\this>lys\tile>lys\tser>lys\tgln>lys\tlys>lys\tmet>lys\tpro>lys\tthr>lys\tphe>lys\tala>lys\tgly>lys\tval>lys\tleu>lys\tasp>lys\targ>lys\ttrp>lys\tglu>lys\ttyr>lys\tcys>met\tasn>met\this>met\tile>met\tser>met\tgln>met\tlys>met\tmet>met\tpro>met\tthr>met\tphe>met\tala>met\tgly>met\tval>met\tleu>met\tasp>met\targ>met\ttrp>met\tglu>met\ttyr>met\tcys>pro\tasn>pro\this>pro\tile>pro\tser>pro\tgln>pro\tlys>pro\tmet>pro\tpro>pro\tthr>pro\tphe>pro\tala>pro\tgly>pro\tval>pro\tleu>pro\tasp>pro\targ>pro\ttrp>pro\tglu>pro\ttyr>pro\tcys>thr\tasn>thr\this>thr\tile>thr\tser>thr\tgln>thr\tlys>thr\tmet>thr\tpro>thr\tthr>thr\tphe>thr\tala>thr\tgly>thr\tval>thr\tleu>thr\tasp>thr\targ>thr\ttrp>thr\tglu>thr\ttyr>thr\tcys>phe\tasn>phe\this>phe\tile>phe\tser>phe\tgln>phe\tlys>phe\tmet>phe\tpro>phe\tthr>phe\tphe>phe\tala>phe\tgly>phe\tval>phe\tleu>phe\tasp>phe\targ>phe\ttrp>phe\tglu>phe\ttyr>phe\tcys>ala\tasn>ala\this>ala\tile>ala\tser>ala\tgln>ala\tlys>ala\tmet>ala\tpro>ala\tthr>ala\tphe>ala\tala>ala\tgly>ala\tval>ala\tleu>ala\tasp>ala\targ>ala\ttrp>ala\tglu>ala\ttyr>ala\tcys>gly\tasn>gly\this>gly\tile>gly\tser>gly\tgln>gly\tlys>gly\tmet>gly\tpro>gly\tthr>gly\tphe>gly\tala>gly\tgly>gly\tval>gly\tleu>gly\tasp>gly\targ>gly\ttrp>gly\tglu>gly\ttyr>gly\tcys>val\tasn>val\this>val\tile>val\tser>val\tgln>val\tlys>val\tmet>val\tpro>val\tthr>val\tphe>val\tala>val\tgly>val\tval>val\tleu>val\tasp>val\targ>val\ttrp>val\tglu>val\ttyr>val\tcys>leu\tasn>leu\this>leu\tile>leu\tser>leu\tgln>leu\tlys>leu\tmet>leu\tpro>leu\tthr>leu\tphe>leu\tala>leu\tgly>leu\tval>leu\tleu>leu\tasp>leu\targ>leu\ttrp>leu\tglu>leu\ttyr>leu\tcys>asp\tasn>asp\this>asp\tile>asp\tser>asp\tgln>asp\tlys>asp\tmet>asp\tpro>asp\tthr>asp\tphe>asp\tala>asp\tgly>asp\tval>asp\tleu>asp\tasp>asp\targ>asp\ttrp>asp\tglu>asp\ttyr>asp\tcys>arg\tasn>arg\this>arg\tile>arg\tser>arg\tgln>arg\tlys>arg\tmet>arg\tpro>arg\tthr>arg\tphe>arg\tala>arg\tgly>arg\tval>arg\tleu>arg\tasp>arg\targ>arg\ttrp>arg\tglu>arg\ttyr>arg\tcys>trp\tasn>trp\this>trp\tile>trp\tser>trp\tgln>trp\tlys>trp\tmet>trp\tpro>trp\tthr>trp\tphe>trp\tala>trp\tgly>trp\tval>trp\tleu>trp\tasp>trp\targ>trp\ttrp>trp\tglu>trp\ttyr>trp\tcys>glu\tasn>glu\this>glu\tile>glu\tser>glu\tgln>glu\tlys>glu\tmet>glu\tpro>glu\tthr>glu\tphe>glu\tala>glu\tgly>glu\tval>glu\tleu>glu\tasp>glu\targ>glu\ttrp>glu\tglu>glu\ttyr>glu\tcys>tyr\tasn>tyr\this>tyr\tile>tyr\tser>tyr\tgln>tyr\tlys>tyr\tmet>tyr\tpro>tyr\tthr>tyr\tphe>tyr\tala>tyr\tgly>tyr\tval>tyr\tleu>tyr\tasp>tyr\targ>tyr\ttrp>tyr\tglu>tyr\ttyr>tyr\n"
column_index_aatypes_transitions = "Species\taromatics>aromatics\tpolar>aromatics\tunpolar>aromatics\tcharged>aromatics\taromatics>polar\tpolar>polar\tunpolar>polar\tcharged>polar\taromatics>unpolar\tpolar>unpolar\tunpolar>unpolar\tcharged>unpolar\taromatics>charged\tpolar>charged\tunpolar>charged\tcharged>charged\n"

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

for p in pairlist: #pairs analysis
  concat.seek(0)
  while 1:
    line=concat.readline()
    seq=concat.readline()
    if p[0] in line:
      seq1=seq[:-1].lower()
    elif p[1] in line:
      seq2=seq[:-1].lower()
    if not line:
      break
  
  if len(p)>2: #bootstrap simulations if the background changed
    length=p[2]
    iterration=p[3]
    plusminus=p[4]
    species=p[5]
    applause='background: '+str(length)+' '+str(iterration)+' '+str(plusminus)+' '+' '.join(species)
    print applause
    stringcounts.append(applause+'\n\n\n')
    stringbiases.append(applause+'\n\n\n')
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

    codons_counts.write(p[0])
    codons_counts.write("\t")
    for value in codonscount.values():
      codons_counts.write(str(value) + "\t")
    codons_counts.write("\t\n")

    codons_counts.write(p[0])
    codons_counts.write("_pvalue\t")
    for value in codonscountpvalue.values():
      codons_counts.write(str(value) + "\t")
    codons_counts.write("\t\n")

    aa_counts.write(p[0])
    aa_counts.write("\t")
    for value in aacount.values():
      aa_counts.write(str(value) + "\t")
    aa_counts.write("\t\n")

    aa_counts.write(p[0])
    aa_counts.write("_pvalue\t")
    for value in aacountpvalue.values():
      aa_counts.write(str(value) + "\t")
    aa_counts.write("\t\n")

    aatypes_counts.write(p[0])
    aatypes_counts.write("\t")
    for value in aaclassifcount.values():
      aatypes_counts.write(str(value) + "\t")
    aatypes_counts.write("\t\n")

    aatypes_counts.write(p[0])
    aatypes_counts.write("_pvalue\t")
    for value in aaclassifcountpvalue.values():
      aatypes_counts.write(str(value) + "\t")
    aatypes_counts.write("\t\n")

    gc_counts.write(p[0])
    gc_counts.write("\t")
    gc_counts.write(str(GC3)+"\t"+str(GC12)+"\t"+str(IVYWREL)+"\t"+str(EKQH)+"\t"+str(PAYRESDGM)+"\t"+str(purineload)+"\t"+str(CvP))
    gc_counts.write("\t\n")

    # end writing
    
    stringcounts=stringcounts[:-1]+[''.join([stringcounts[-1],("counting of %s\n\n" % p[0])+''.join(strcountings(codonscount, aacount, aaclassifcount,codonscountpvalue,aacountpvalue,aaclassifcountpvalue,GC3,GC12,IVYWREL,EKQH,PAYRESDGM,purineload,CvP))+'\n\n'])]
  else:
    towrite=1
  # -----------------------------

  stringbiases.append("mutation biases from %s to %s\n\n" % (p[1], p[0]))
  stringbiases=stringbiases+strbiases(codons, aa, aaclassif,codonspvalue,aapvalue,aaclassifpvalue,aa_transitions,aatypes_transitions,p)
  stringbiases.append('\n\n') 
  print 'done'

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
