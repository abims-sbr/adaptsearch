#!/usr/bin/env python
#keeps the longest ORF found in the 6 possible ORF alltogether
#python find_ORF.py file output

def find_orf(entry):
  orf={}
  orf_length={}
  stop=['TAA','TAG','TGA']
  for i in range(0,3):
    pos=i
    orf[i]=[0]
    while pos<len(entry):
      if entry[pos:pos+3] in stop:
        orf[i].append(pos-1)
        orf[i].append(pos+3)
      pos+=3
    orf[i].append(len(entry)-1)
    orf_length[i]=[]
    for u in range(1,len(orf[i])):
      orf_length[i].append(orf[i][u]-orf[i][u-1]+1)
    orf[i]=[orf[i][orf_length[i].index(max(orf_length[i]))],orf[i][orf_length[i].index(max(orf_length[i]))+1]]
  orf_max={0:max(orf_length[0]),1:max(orf_length[1]),2:max(orf_length[2])}
  orf=orf[max(list(orf_max.keys()), key=(lambda k: orf_max[k]))]
  if orf[0]==0:
    orf[0]=orf[0]+max(list(orf_max.keys()), key=(lambda k: orf_max[k]))
  return orf


def reverse_seq(entry):
  nt={'A':'T','T':'A','G':'C','C':'G', 'N':'N'}
  seqlist=[]
  for i in range(len(entry)-1,-1,-1):
    seqlist.append(nt[entry[i]])
  seq=''.join(seqlist)
  return seq
  


# RUN

import string, os, sys, re

path_IN = sys.argv[1]
file_OUT = open(sys.argv[2], "w")
f_in = open(path_IN, "r")
inc=1
threshold=0 #minimal length of the ORF
while 1:
  line=f_in.readline()
  if not line:
    break
  line=f_in.readline()
  name=">"+path_IN[:2]+str(inc)+"_1/1_1.000_"
  high_plus=find_orf(line[:-1])
  reverse=reverse_seq(line[:-1])
  high_minus=find_orf(reverse)
  if high_plus[1]-high_plus[0]>threshold or high_minus[1]-high_minus[0]>threshold:
    inc+=1
    if high_plus[1]-high_plus[0]>high_minus[1]-high_minus[0]:
      file_OUT.write("%s" %name)
      file_OUT.write(str(high_plus[1]-high_plus[0]+1)+"\n")
      file_OUT.write("%s" %line[high_plus[0]:high_plus[1]+1])
      file_OUT.write("\n")
    else:
      file_OUT.write("%s" %name)
      file_OUT.write(str(high_minus[1]-high_minus[0]+1)+"\n")
      file_OUT.write("%s" %reverse[high_minus[0]:high_minus[1]+1])
      file_OUT.write("\n")
f_in.close()
file_OUT.close()

  
