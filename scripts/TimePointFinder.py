"""
@author: hocinerekaik

Find corresponding time-point of a given dataset when compared to computed timepoint matrix generated with BigWigMatrixGen.py

Modified May 20 2022 by lldelisle
To propose command line and dynamic norm chrom setting
"""

import numpy as np
import pyBigWig
import re
import tkinter as tk
from tkinter import filedialog, simpledialog
import argparse
import os

argp = argparse.ArgumentParser(
    description=("Find corresponding time-point of a given dataset when compared to computed timepoint matrix generated with BigWigMatrixGen.py"))
argp.add_argument('--normChr', default=None,
                  help="which chromosome should be used to evaluate the time point (chr6=HoxA).")
argp.add_argument('--bed', default=None,
                  help="bed file with regions on which coverage is computed")
argp.add_argument('--bigwig', default=None,
                  help="bigwig file for which the time point will be computed")
argp.add_argument('--matrix', default=None,
                  help="matrix file obtained by BigWigMatrixGen.py")
args = argp.parse_args()

# Initialize root
root = None
lFilesList = []

# cluster used for normalization chr6 = HoxA
normChr = args.normChr
if normChr is None:
    # If no argument is given use interactive mode
    root = tk.Tk()
    root.withdraw()
    # Let the user choose the chromosome
    normChr = simpledialog.askstring("Chromosome norm selection",
                                     'Enter the chromosome you want to use to find the time point:' +
                                     '\n chr6 for HoxA')

"""
Bed file
"""

# bed file intervals covering regions of interest. The same used in BigWigMatrixGen.py
if args.bed is None:
    if root is None:
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

with open(sFPath_Bed) as f:
    lBinsLines = f.read().splitlines()


"""
Open bigwig file and get coverage
"""

if args.bigwig is None:
    # If no argument is given use interactive mode
    """
    Select bigwig file
    """
    if root is None:
        root = tk.Tk()
        root.withdraw()

    sFPath_Bigwig = filedialog.askopenfilename(title="Select BigWig file ")
    root.update()
    if sFPath_Bigwig == "":
        exit("Please provide a BigWig file")
else:
    sFPath_Bigwig = args.bigwig
    if not os.path.exists(sFPath_Bigwig):
        raise Exception(sFPath_Bigwig + " does not exist")

bw = pyBigWig.open(sFPath_Bigwig)
lRawTimePoint = []
for line in lBinsLines:
    lSplittedLine = line.split("\t")
    # Only compute coverage in the normChr
    if lSplittedLine[0] == normChr:
        lRawTimePoint.append([lSplittedLine[0], lSplittedLine[1], lSplittedLine[2], bw.stats(
            lSplittedLine[0], int(lSplittedLine[1]), int(lSplittedLine[2]), type="mean")[0]])

bw.close()

"""
Compare and find the best timepoint
"""
# Get the matrix:

if args.matrix is None:
    # If no argument is given use interactive mode
    """
    Select matrix file
    """
    if root is None:
        root = tk.Tk()
        root.withdraw()

    sFPath_Matrix = filedialog.askopenfilename(title="Select Matrix file ")
    root.update()
    if sFPath_Bigwig == "":
        exit("Please provide a Matrix file")
else:
    sFPath_Matrix = args.matrix
    if not os.path.exists(sFPath_Matrix):
        raise Exception(sFPath_Matrix + " does not exist")

with open(sFPath_Matrix) as fMatrix:
    lMatrixLines = fMatrix.read().splitlines()

# Subset lMatrixLines to get only the lines corresponding to the chromosome
lMatrix = []
for line in lMatrixLines:
    lSplittedLine = line.split("\t")
    if lSplittedLine[0] == normChr:
        lMatrix.append([lSplittedLine[0], lSplittedLine[1],
                        lSplittedLine[2], lSplittedLine[3], float(lSplittedLine[4])])

# We asume the matrix has exactly the same intervals on the norm.chr
# We asume the matrix is:
# - for each time point
#   - for each interval of the bed file
#      - the coverage
# lMatrix[k + len(lRawTimePoint) * i][4] is the coverage
# for interval k at time point of indice i

# Initialize minDelta with no time and a very high value
minDelta = ['', 10000000000.0]
for i in range(0, int(len(lMatrix) / len(lRawTimePoint)), 1):
    # i will be an index for a time point in the lMatrix
    delta = 0.0
    working_time_point = lMatrix[len(lRawTimePoint)*i][3]
    for NumberInterval, interval in enumerate(lRawTimePoint):
        # inverval is like [chr, start, end, average]
        index_in_lMatrix = NumberInterval + (len(lRawTimePoint) * i)
        # Check this is the same interval
        assert interval[0] == lMatrix[index_in_lMatrix][0]
        assert interval[1] == lMatrix[index_in_lMatrix][1]
        assert interval[2] == lMatrix[index_in_lMatrix][2]
        # And the same time point
        assert lMatrix[index_in_lMatrix][3] == working_time_point
        if interval[3] == None:
            interval[3] = 0
        delta = delta + \
            abs(interval[3] - lMatrix[index_in_lMatrix][4])
    if minDelta[1] > delta:
        minDelta = [working_time_point, delta]
print(minDelta[0] + " " + str(minDelta[1]))
