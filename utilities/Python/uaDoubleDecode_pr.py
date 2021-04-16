pass
# uaDoubleDecode_pr.py
#
# This doubly decodes User-Agent strings in data files from Philipp Richter,
# and writes them back out again.
# 
# Python 3, run in Enthought Python (https://www.enthought.com/).
# Based upon version for Anaconda Python first written 2016-2017.
#
# Last changed 16th August 2018.

import gzip
import os
import sys
import numpy as np
import scipy as sp
#import pandas as pd
import urllib.parse
import string as st
import fnmatch
import pickle as pik
import functools as ft
import random as ra
import math as ma
import itertools as it
import unicodedata as ud

#in_path = "c:/builds/2017/ua/ua-1h/"
in_path = "c:/builds/2018/UA/PR/"
out_path = "c:/builds/2018/UA/decoded/"

rejectsFile = "c:/builds/2018/UA/rejects/pr_rejects.data"

#included_pattern = '*_pfx.gz'
included_pattern = '*.gz'

rejectOut = open(rejectsFile, 'w', encoding='utf8')

FileNames = []
for dirpath, dirnames, filenames in os.walk(in_path):
    for filename in fnmatch.filter(filenames, included_pattern):
        FileNames.append(os.path.join(dirpath, filename))
# From http://stackoverflow.com/questions/16229982/read-large-gzip-files-in-python

fileCount = 0
period = 2000000
totalLineCount = 0
linesWritten = 0
linesRejected = 0
linesSkipped = 0

for fName in FileNames:
    filename_w_ext = os.path.basename(fName)
    filename, file_extension = os.path.splitext(filename_w_ext)
    with gzip.open(fName, 'r') as infile:
        outFile1 = "%s%s.data"%(out_path, filename)
        pOut1 = open(outFile1, 'w', encoding='utf8')
#   with open(fName, 'r') as infile:
        lineCountWithinFile = 0
        print("Doing file '%s' ..."%fName)
        fileCount += 1
        while True:
            try:
                line = infile.readline()
            except:
                break
            lineStr = str(line, encoding='utf8')
#           lineStr = line # No UTF-8 decoding needed
            lineStr = lineStr.rstrip('\n')
            if 0 == len(lineStr):
                break
            totalLineCount += 1
            lineCountWithinFile += 1
            reject = False
            reason = ""
            skipped = False
            if lineStr.startswith('#'):
                skipped = True
            else:
                fields = lineStr.split(sep='\t')
                if 7 == len(fields):
                    extracted = fields[6]
                    rewrite = urllib.parse.unquote(urllib.parse.unquote(extracted))
                    rewrite = rewrite.replace("\n", "<<newline>>").replace("\t", "<<tab>>")
                else:
                    reason = "failed match"
                    rewrite = ""
                    reject = True
            if skipped:
                    linesSkipped += 1
            elif reject:
                rejectOut.write("Rejected line %s in file '%s' tc %s, reason '%s', extract: '%s':\n'%s'"%
                        (lineCountWithinFile, fName, totalLineCount, reason, rewrite, lineStr))
                linesRejected += 1
            else:
                try:
                    pOut1.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(fields[0], fields[1], fields[2], fields[3],
                                                                fields[4], fields[5], rewrite))
                    linesWritten += 1
                except:
                    reason = "unsupported coding"
                    rejectOut.write("Rejected line %s in file '%s' tc %s, reason '%s', extract: '%s':\n'%s'"%
                                    (lineCountWithinFile, fName, totalLineCount, reason, rewrite, lineStr))
                    linesRejected += 1
            if 0 == (totalLineCount%period):
                print("... did file %s, '%s', %s lines within file, written %s, %s lines rejected, %s skipped ..."%
                          (fileCount, fName, lineCountWithinFile, linesWritten, linesRejected, linesSkipped))
        infile.close()
        pOut1.close()
        print("... did file %s, '%s', %s lines within file, written %s, %s lines rejected, %s skipped ..."%
                          (fileCount, fName, lineCountWithinFile, linesWritten, linesRejected, linesSkipped))
#       break

rejectOut.close()

print("Completed reading %s files, %s total lines, written %s, %s lines rejected, %s lines skipped.\n"%
      (fileCount, totalLineCount, linesWritten, linesRejected, linesSkipped))



