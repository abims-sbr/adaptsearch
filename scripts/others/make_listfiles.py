#!/usr/bin/env python

# This scripts writes a config file containing all the file names specified in the command line (one file per line)
# Victor Mataigne - 06/18

import sys

lfiles = sys.argv[1].split(',')
list_files = open('list_files', 'w')

for file in lfiles:
    list_files.write(file+'\n')

list_files.close()