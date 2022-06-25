#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script compute the cumulative coverage of given bigwig files at regions of interest. The output is a bedgraph file with the new computed values. 
"""


"""
Created on Thu Apr 16 13:48:43 2020

@author: hocinerekaik

Modified May 9 2022 by lldelisle
To propose command line
"""

import numpy as np
import pyBigWig
import random
import re
import argparse
import os


argp = argparse.ArgumentParser(
    description=("Compute cumulative coverage on defined regions"))
argp.add_argument('--bigwig', default=None, nargs='+',
                  help="bigwig files (separated by space)")
argp.add_argument('--bed', default=None,
                  help="bed file with regions on which coverage is computed")
argp.add_argument('--output', default=None,
                  help="output bedgraph file")

args = argp.parse_args()

if None in [args.bigwig, args.bed]:
    import tkinter as tk
    from tkinter import filedialog

# Initialize root
root = None

"""
Bed file
"""
# bed file with regions of interest

if args.bed is None:
    # If no argument is given use interactive mode
    root = tk.Tk()
    root.withdraw()
    sFPath_Bed = filedialog.askopenfilename(
        title="Select bed file to compute coverage on")
    root.update()
else:
    sFPath_Bed = args.bed
    # Check files exists:
    if not os.path.exists(sFPath_Bed):
        raise Exception(sFPath_Bed +  " does not exist")

with open(sFPath_Bed, 'r') as f:
    lBinsLines = f.read().splitlines()


"""
Get Bigwig files
"""

lFilesList = []
if args.bigwig is None:
    # If no argument is given use interactive mode
    if root is None:
        root = tk.Tk()
        root.withdraw()
    sFPath_Bigwig = "NA"
    i = 1
    while sFPath_Bigwig != "":
        sFPath_Bigwig = filedialog.askopenfilename(
            title="Select Bigwig file " + str(i))
        root.update()
        if sFPath_Bigwig == "":
            break
        lFilesList.append(sFPath_Bigwig)
        i = i + 1
else:
    # Check files exists:
    for bw in args.bigwig:
        if not os.path.exists(bw):
            raise Exception(bw + " does not exist")
        lFilesList.append(bw)

"""
Open bigwig file and get coverage
"""
mBigWigsCoverage = []
# Will be an array with one element per bigwig
# Each element is an array with one element per line in bed file
# Each element is [chr, start, end, average]
for sample in lFilesList:
    # sample is a file name of bigwig
    bw = pyBigWig.open(sample)
    lRawTimePoint = []
    for line in lBinsLines:
        lSplittedLine = line.split("\t")
        lRawTimePoint.append([lSplittedLine[0], lSplittedLine[1], lSplittedLine[2], bw.stats(
            lSplittedLine[0], int(lSplittedLine[1]), int(lSplittedLine[2]), type="mean")[0]])
    mBigWigsCoverage.append(lRawTimePoint)


"""
Compute cumutative Coverage
"""
ModifiedIntervals = []
# Will be an array with one element per line in bed file
# Each element is [chr, start, end, sum of average]
for x, line in enumerate(lBinsLines):
    lSplittedLine = line.split("\t")
    ModifiedIntervals.append(
        [lSplittedLine[0], lSplittedLine[1], lSplittedLine[2], 0.0])
    for sample in mBigWigsCoverage:
        if sample[x][3] == None:
            sample[x][3] = 0.0
        ModifiedIntervals[x][3] = ModifiedIntervals[x][3] + sample[x][3]


"""
save results in bedgraph file
"""
if args.output is None:
    output_fn = lFilesList[0][0:len(sFPath_Bigwig) - 7] + '_CumBW.bedgraph'
else:
    output_fn = args.output
fResult = open(output_fn, 'w')

for Interval in ModifiedIntervals:
    fResult.write("%s\n" % (Interval[0] + "\t" +
                            Interval[1] + "\t" + Interval[2] + "\t" + str(Interval[3])))

fResult.close()
