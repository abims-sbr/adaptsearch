#!/usr/bin/env python

import sys

lfiles = sys.argv[1].split(',')
list_files = open('list_files', 'w')

for file in lfiles:
    list_files.write(file+'\n')

list_files.close()