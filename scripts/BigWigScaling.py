#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 13:48:43 2020

@author: hocinerekaik

Modified May 4 2022 by lldelisle
To propose command line

"""

import numpy as np
import pyBigWig
import random
import argparse
import os

argp = argparse.ArgumentParser(
    description=("Normalize bigwigs using corresponding peaks and min max method"))
argp.add_argument('--referenceBW', default=None,
                  help="reference bigwig")
argp.add_argument('--referenceNP', default=None,
                  help="reference narrowPeak")
argp.add_argument('--otherBW', default=None, nargs='*',
                  help="other bigwigs to normalize separated by space.")
argp.add_argument('--otherNP', default=None, nargs='*',
                  help="other narrowPeaks corresponding to otherBW (same order).")
argp.add_argument('--seed', default=None, type=int,
                  help="set seed for reproducibility")
argp.add_argument('--prefixReports', default="",
                  help="prefix to use for reports")
argp.add_argument('--noOutput', action="store_true",
                  help="For debug purpose, do not output the normalized bigwig")
args = argp.parse_args()

if None in [args.referenceBW, args.referenceNP, args.otherBW, args.otherNP]:
    import tkinter as tk
    from tkinter import filedialog

# Generate a list with macs2 file path, bigwig file path,
# min = average in interpeak regions and
# max = average of max coverage in common peaks
lFilesList = []

# If no argument were given use interactive mode
if args.referenceBW is None:
    """
    Open peaks and Bigwig files
    """
    root = tk.Tk()
    root.withdraw()


    sFPath_Macs2 = "NA"
    i = 1
    while sFPath_Macs2 != "":
        sFPath_Macs2 = filedialog.askopenfilename(
            title="Select Macs2 peaks file " + str(i))
        root.update()
        if sFPath_Macs2 == "":
            break
        sFPath_Bigwig = filedialog.askopenfilename(
            title="Select Bigwig file " + str(i))
        root.update()
        lFilesList.append([sFPath_Macs2, sFPath_Bigwig, 0.0, 0.0])
        i = i + 1
else:
    # Perform some checks
    if not os.path.exists(args.referenceNP):
        raise Exception(args.referenceNP + " does not exists.")
    if not os.path.exists(args.referenceBW):
        raise Exception(args.referenceBW + " does not exists.")
    if args.otherBW is None:
        raise Exception("No otherBW specified.")
    if args.otherNP is None:
        raise Exception("No otherNP specified.")
    if len(args.otherBW) != len(args.otherNP):
        raise Exception("otherBW and otherNP does not contain the same number of items.")
    # Start with reference:
    lFilesList.append([args.referenceNP, args.referenceBW, 0.0, 0.0])
    # Continue with other:
    for peak, bw in zip(args.otherNP, args.otherBW):
        lFilesList.append([peak, bw, 0.0, 0.0])

# Load in memory all control peaks
with open(lFilesList[0][0], "r") as fpeakcontrol:
    print("MACS 2 peaks control file :", lFilesList[0][0])
    lPeaksControlLines = fpeakcontrol.read().splitlines()
    lastline = len(lPeaksControlLines)


"""
Open all files in lFilesList
"""
lOpenedFiles = []
lMaxPeaksFile = []
MaxPeaksHeader = lFilesList[0][0]
for sample in lFilesList:
    fpeaksample = open(sample[0], "r")
    bw = pyBigWig.open(sample[1])
    lOpenedFiles.append([fpeaksample, bw])
    MaxPeaksHeader = MaxPeaksHeader + "\t" + str(sample[0])


"""
Verify if the random line in fpeakcontrol is present in all compared files. If Macs2 reported the peak in all samples.
"""
# Change number of common peaks needed
numberTest = 1000
lMaxPeaksFile = [None] * (numberTest + 1)
lMaxPeaksFile[0] = MaxPeaksHeader
i = 0
lValidatedLines = []
# If specified, use seed
if args.seed is not None:
    random.seed(args.seed)
    np.random.seed(args.seed)

# i is the index of the current peak to find
while i != numberTest:
    # x is the line number in the control peak
    x = random.randint(0, lastline - 1)
    ControlPeakLine = lPeaksControlLines[x]
    lSplittedControlLine = ControlPeakLine.split("\t")
    s = 1
    for sample in lOpenedFiles:
        # sample is [narrowPeak_opened_file, bigwig_opened]
        # Come back to the beginning of the file
        sample[0].seek(0)
        # Get all lines
        lSamplePeaksLines = sample[0].readlines()
        # Store all compatible peaks
        lCandidateSamplePeaks = []
        # For each line in peak
        # Check if start is less than 5kb from the control peak line
        for y in range(0, len(lSamplePeaksLines)):
            lSplittedSampleLine = lSamplePeaksLines[y].split("\t")
            if lSplittedSampleLine[0] == lSplittedControlLine[0] and abs(int(lSplittedSampleLine[1]) - int(lSplittedControlLine[1])) < 5000:
                lCandidateSamplePeaks.append(lSplittedSampleLine)
        # If more than one peaks are found
        if len(lCandidateSamplePeaks) > 0:
            # Only keep the closest from the start
            BestCandidate = lCandidateSamplePeaks[0]
            for Candidate in lCandidateSamplePeaks:
                if abs(int(Candidate[1]) - int(lSplittedControlLine[1])) < abs(int(BestCandidate[1]) - int(lSplittedControlLine[1])):
                    BestCandidate = Candidate
            # If the best candidate has start and end closer than 400bp of start and end of the control peak line
            # We consider that it is valid and store it in lMaxPeaksFile
            if abs(int(BestCandidate[1]) - int(lSplittedControlLine[1])) < 400 and abs(int(BestCandidate[2]) - int(lSplittedControlLine[2])) < 400:
                CandidateFound = True
                if s == 1:
                    lMaxPeaksFile[i+1] = "\t" + BestCandidate[0] + \
                        ":" + BestCandidate[1] + "-" + BestCandidate[2]
                else:
                    lMaxPeaksFile[i+1] = str(lMaxPeaksFile[i+1]) + "\t" + BestCandidate[0] + \
                        ":" + BestCandidate[1] + "-" + BestCandidate[2]
                s = s + 1
            else:
                CandidateFound = False
        else:
            CandidateFound = False
    # If a candidate was found in all samples:
    if CandidateFound == True:
        # Increment i (index of peak to find)
        i = i + 1
        # Store the line number x
        lValidatedLines.append(x)
print(str(numberTest) + " homologous peaks found")

# Close all files:
for sample in lOpenedFiles:
    # sample is [narrowPeak_opened_file, bigwig_opened]
    sample[0].close()
    sample[1].close()

"""
Start batch treatment of all samples
"""
for n, sample in enumerate(lFilesList):
    """
    Initiate BigWig and Macs2 files 
    """
    fpeaksample = open(sample[0], "r")
    bw = pyBigWig.open(sample[1])
    if sample[1][len(sample[1]) - 6:len(sample[1])] == "bigwig":
        sFPath_NewBigwig = sample[1][0:len(
            sample[1]) - 7] + "_Normalized.bigwig"
    else:
        sFPath_NewBigwig = sample[1] + "_Normalized.bigwig"
    print("BigWig file :", sample[1])
    if not args.noOutput:
        bwNew = pyBigWig.open(sFPath_NewBigwig, "w")
        lNewheader = d = [(k, v) for k, v in bw.chroms().items()]
        bwNew.addHeader(lNewheader, maxZooms=4)

    lPeaksLines = fpeaksample.read().splitlines()
    lastline = len(lPeaksLines)

    """
    Find mean coverage of inter-peaks regions
    """
    # Generate random numberTest x 2000bp intervals in inter-peaks regions
    # First, generate numpy array with random lines, to be applied to peaks file
    numberTest = 100
    aMeanInterval = np.empty(numberTest)
    aMeanInterval[:] = np.NaN
    # Take randomly numberTest indices
    # Corresponding to peaks
    randomLines = np.random.rand(numberTest, 1) * lastline
    randomLines = randomLines.astype(int)
    # Retrieve intervals from inter-peaks regions
    for x, i in enumerate(np.nditer(randomLines)):
        # x is the index of the random interval
        # i is the index of the random line from which we look for an interpeak region
        endFirstPeak = 0
        startSecondPeak = 0
        j = 0
        chr_ipj = "chr"
        chr_ipjp1 = "chr"
        # j is the offset from i and will increase until there is
        # more than 10kb between peak (i+j) and peak (i+j+1)
        # Verify if distance bewteen peak1 and peak2 is > 10000 bp : avoid to catch intervals with rich coverage
        # Verify that consecutive peaks are still on the same chromosome
        # And chromosome name does not contains '_' (contig)
        while chr_ipj != chr_ipjp1 or '_' in chr_ipj or (startSecondPeak - endFirstPeak) < 100000:
            # If end of file reached of consecutive peaks are on different chromosomes
            if i+j+1 >= lastline or chr_ipj != chr_ipjp1 or '_' in chr_ipj:
                # Take another random value
                i = random.randint(0, lastline - 1)
                j = 0
            strFirstPeakLine = lPeaksLines[i+j].split("\t")
            endFirstPeak = int(strFirstPeakLine[2])
            strSecondPeakLine = lPeaksLines[i+j+1].split("\t")
            startSecondPeak = int(strSecondPeakLine[1])
            chr_ipj = strFirstPeakLine[0]
            chr_ipjp1 = strSecondPeakLine[0]
            j = j + 1
        # Compute coverage in this interval:
        # From the middle of the interpeak to 2kb on the right.
        startInterval = int(((startSecondPeak-endFirstPeak)/2)+endFirstPeak)
        aMeanInterval[x] = bw.stats(
            strFirstPeakLine[0], startInterval, startInterval + 2000)[0]
    # Compute mean coverage accross these inter-peak regions
    print(np.nanmean(aMeanInterval))
    minMeanInterval = np.nanmean(aMeanInterval, dtype='float64')
    lFilesList[n][2] = minMeanInterval

    """
    Find max coverage of MACS2 peaks
    """""
    # Generage random x numberTest intervals from MACS2 peaks list computed previously
    # Retrieve max coverage of each peak and compute mean coverage of those maximum values
    x = 0
    npMaxPeaks = np.empty(len(lValidatedLines))
    for i in lValidatedLines:
        strPeakLine = lPeaksControlLines[i].split("\t")
        npMaxPeaks[x] = bw.stats(strPeakLine[0], int(
            strPeakLine[1]), int(strPeakLine[2]), type="max")[0]
        x = x + 1
    print(np.nanmean(npMaxPeaks))
    maxMeanPeaks = np.nanmean(npMaxPeaks)
    lFilesList[n][3] = maxMeanPeaks

    """
    Core function : For each base, compute new value usin the Rescaling method. x' = (x - min(x)) / (max(x) - min(x))
    """

    def fScaling():
        difPeaksIntervals = (maxMeanPeaks - minMeanInterval) / 100
        for i in bw.chroms():
            # For each chromosome
            print(i, "..Scaling")
            for interval in bw.intervals(i):
                # For each entry in the original bigwig
                # Keep same start end
                start = interval[0]
                end = interval[1]
                # Use normalized value
                NewValue = (interval[2] - minMeanInterval) / difPeaksIntervals
                # Avoid negative value
                if NewValue < 0:
                    NewValue = 0.0
                bwNew.addEntries([i], [start], ends=[end], values=[NewValue])
            print(i, "..Processed")

    if not args.noOutput:
        # Run the scaling
        fScaling()
        # Close the new bigwig
        bwNew.close()
    # Close the original bigwig
    bw.close()
    print("Job finished")

# Write statistics file
dirname = os.path.dirname(args.prefixReports)
if dirname != '' and not os.path.isdir(dirname):
    print(dirname + " does not exists. Will write to current directory.")
    args.prefixReports = os.path.basename(args.prefixReports)

with open(args.prefixReports + 'Stats.txt', 'w') as fResult:
    fResult.write("%s\n" % ("MACS2 peaks file" + "\t" +
                            "BigWig file" + "\t" + "Min value" + "\t" + "Max value"))
    for item in lFilesList:
        fResult.write("%s\n" % (str(item[0]) + "\t" +
                                str(item[1]) + "\t" + str(item[2]) + "\t" + str(item[3])))


with open(args.prefixReports + 'Max peaks.txt', 'w') as fResult:
    for item in lMaxPeaksFile:
        fResult.write("%s\n" % item)
