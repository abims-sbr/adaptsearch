#!/usr/bin/python
#launched by correctcodeml.sh

import string, os, sys, re

path_IN=sys.argv[1]
#mgene=sys.argv[2]
mgene=0
folder=string.split(sys.argv[3],'/')[::-1][1]
values=string.split(path_IN,"/")
values.pop()
values=values[::-1]
output=[]
nsvalue=0
ff='free'
model=0

i=0
while folder not in values[i]:
  output.append(values[i])
  if 'ns' in values[i]:
    nsvalue=values[i][2]
  if 'w' in values[i]:
    ff=values[i][1:]
  if 'model' in values[i]:
    model=values[i][-1:]
  i+=1
    
outfile='_'.join(output[::-1])

files=os.listdir(path_IN)
for i in files:
  if i[-3:]=='txt':
    tree=i
  if (i[-3:]=='phy') or (i[-5:]=='fasta'):
    seq=i


codeml=open(path_IN+"codeml.ctl","w")
codeml.write("      seqfile = %s * sequence data file name \n" % seq )
codeml.write("      outfile = %s * main result file name \n" % outfile)
codeml.write("     treefile = %s * tree structure file name\n" % tree)
codeml.write("        noisy = 9  * 0,1,2,3,9: how much rubbish on the screen \n      verbose = 0  * 1: detailed output, 0: concise output \n      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic \n                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise \n      seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs \n    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table \n*       ndata = 10 \n        clock = 0   * 0:no clock, 1:clock; 2:local clock \n       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a \n                   * 7:AAClasses \n   aaRatefile = wag.dat * only used for aa seqs with model=empirical(_F) \n                  * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own \n")
codeml.write("        model = %s \n                   * models for codons: \n                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches \n                   * models for AAs or codon-translated AAs: \n                       * 0:poisson, 1:proportional,2:Empirical,3:Empirical+F \n                       * 6:FromCodon, 8:REVaa_0, 9:REVaa(nr=189) \n" %model)
codeml.write("      NSsites = %s  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs; \n                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma; \n                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1; \n                   * 13:3normal>0 \n        icode = 0  * 0:universal code; 1:mammalian mt; 2-11:see below \n        Mgene = %s  * 0:rates, 1:separate;  \n    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated \n        kappa = 2  * initial or fixed kappa\n" % (nsvalue,mgene))
if ff=="fixed":
  codeml.write("    fix_omega = 1  * 1: omega or omega_1 fixed, 0: estimate  \n        omega = 1.0 * initial or fixed omega, for codons or codon-based AAs \n    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha \n        alpha = 0. * initial or fixed alpha, 0:infinity (constant rate) \n       Malpha = 0  * different alphas for genes \n        ncatG = 3  * # of categories in dG of NSsites models \n      fix_rho = 1  * 0: estimate rho; 1: fix it at rho \n          rho = 0. * initial or fixed rho,   0:no correlation \n        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates \n RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2) \n   Small_Diff = .5e-6 \n*   cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)? \n* fix_blength = 0   * 0: ignore, -1: random, 1: initial, 2: fixed \n       method = 0   * 0: simultaneous; 1: one branch at a time ")
if ff=="free":
  codeml.write("    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate  \n        omega = 0.2 * initial or fixed omega, for codons or codon-based AAs \n    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha \n        alpha = 0. * initial or fixed alpha, 0:infinity (constant rate) \n       Malpha = 0  * different alphas for genes \n        ncatG = 3  * # of categories in dG of NSsites models \n      fix_rho = 1  * 0: estimate rho; 1: fix it at rho \n          rho = 0. * initial or fixed rho,   0:no correlation \n        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates \n RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2) \n   Small_Diff = .5e-6 \n*   cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)? \n* fix_blength = 0   * 0: ignore, -1: random, 1: initial, 2: fixed \n       method = 0   * 0: simultaneous; 1: one branch at a time ")

codeml.close()
