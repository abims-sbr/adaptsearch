#!/usr/bin/env python

import zipfile, os, re

f_25 = "^25_.*$"
f_19 = "^19_Reciprocal.*$"

f_DNA = zipfile.ZipFile("output_file_DNA.zip", "w")
f_PROT = zipfile.ZipFile("output_file_PROT.zip", "w")

a = os.listdir("./")
a = sorted(a)

for i in a :
	if re.match(f_25, i) :
		f_DNA.write(i)
	if re.match(f_19, i) :
		f_PROT.write(i)
