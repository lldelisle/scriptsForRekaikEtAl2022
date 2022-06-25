"""
@author: hocinerekaik

This script computes the variation (extracted from the matrice produced with BigWigMatrixGen.py) between a target time-point and the computed time-point (found with TimePointFinder.py) for each bin,
and applies the changes on the input bigwig file. The output is a bedgraph with the modifed coverage value for all clusters.

Modified May 20 2022 by lldelisle
To propose command line
"""

import numpy as np
import pyBigWig
import re
import tkinter as tk
from tkinter import filedialog
import argparse
import os

argp = argparse.ArgumentParser(
    description=("This script computes the variation"
                 " (extracted from the matrice produced "
                 "with BigWigMatrixGen.py) between a target "
                 "time-point and the computed time-point"
                 " (found with TimePointFinder.py) for each bin,"
                 "and applies the changes on the input bigwig file."
                 " The output is a bedgraph with the modifed coverage"
                 " value for all clusters."))
argp.add_argument('--computedTimePoint', default=None, type=int,
                  help="Time point found by TimePointFinder.py")
argp.add_argument('--targetTimePoint', default=None, type=int,
                  help="Time point objective")
argp.add_argument('--bed', default=None,
                  help="bed file with regions on which coverage is computed")
argp.add_argument('--bigwig', default=None,
                  help="bigwig file to normalize")
argp.add_argument('--matrix', default=None,
                  help="matrix file obtained by BigWigMatrixGen.py")
argp.add_argument('--normWindowSize', default=None, type=int,
                  help="Number of intervals to consider to evaluate" +
                  " the normalizing factor (variation).")
args = argp.parse_args()

# Initialize root
root = None
lFilesList = []

ComputedTimePoint = args.computedTimePoint
if ComputedTimePoint is None:
    # If no argument is given use interactive mode
    root = tk.Tk()
    root.withdraw()
    ComputedTimePoint = simpledialog.askinteger("Computed Time point selection",
                                                'Enter the time point found by TimePointFinder.py in hours')


TargetTimePoint = args.targetTimePoint
if TargetTimePoint is None:
    # If no argument is given use interactive mode
    if root is None:
        root = tk.Tk()
        root.withdraw()
    TargetTimePoint = simpledialog.askinteger("Targetted Time point selection",
                                              'Enter the time point to which the bigwig should be scaled in hours')

"""
Bed file
"""
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

# bed file intervals covering regions of interest. The same used in BigWigMatrixGen.py
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
    lRawTimePoint.append([lSplittedLine[0], lSplittedLine[1], lSplittedLine[2], bw.stats(
        lSplittedLine[0], int(lSplittedLine[1]), int(lSplittedLine[2]), type="mean")[0]])

bw.close()

"""
Open the matrix
"""

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

# Store all values into array
lMatrix = []
for line in lMatrixLines:
    lSplittedLine = line.split("\t")
    lMatrix.append([lSplittedLine[0], lSplittedLine[1],
                    lSplittedLine[2], lSplittedLine[3], float(lSplittedLine[4])])


"""
Apply variation
"""

# We asume the matrix has exactly the same intervals on the norm.chr
# We asume the matrix is:
# - for each time point
#   - for each interval of the bed file
#      - the coverage
# lMatrix[k + len(lRawTimePoint) * i][4] is the coverage
# for interval k at time point of indice i

# First find the index corresponding to the TargetTimePoint
IndexTargetTimePoint = 0
for line in lMatrix:
    if line[3] == str(TargetTimePoint):
        break
    IndexTargetTimePoint = IndexTargetTimePoint + 1
# Then find the index corresponding to the ComputedTimePoint
IndexComputedTimePoint = 0
for line in lMatrix:
    if line[3] == str(ComputedTimePoint):
        break
    IndexComputedTimePoint = IndexComputedTimePoint + 1



ModifiedIntervals = []
# The objective is to apply to each interval a correction
# new value = old value / 
#   (value in ref at computed time point / value in ref at target time point)
# The issue is that if we consider only one bin
# The normalizing factor = ref_computed / ref_target
# Is really instable
# Therefore we compute the normalizing factor on a sliding window

# Get the window size for normalizing factor
normWindowSize = args.normWindowSize
if normWindowSize is None:
    # If no argument is given use interactive mode
    if root is None:
        root = tk.Tk()
        root.withdraw()
    normWindowSize = simpledialog.askinteger("Normalization window size selection",
                                             'Enter the number of intervals to consider' +
                                             ' to evaluate the normalizing factor (variation)')
assert normWindowSize > 0, \
    "normWindowSize should be at least 1"

for i, line in enumerate(lRawTimePoint):
    # line is like [chr, start, end, average]
    # for the bigwig to normalize
    target_i_index = IndexTargetTimePoint + i
    computed_i_index = IndexComputedTimePoint + i
    # Check the intervals and time correspond:
    assert line[0] == lMatrix[target_i_index][0]
    assert line[1] == lMatrix[target_i_index][1]
    assert line[2] == lMatrix[target_i_index][2]
    assert str(TargetTimePoint) == lMatrix[target_i_index][3]
    assert line[0] == lMatrix[computed_i_index][0]
    assert line[1] == lMatrix[computed_i_index][1]
    assert line[2] == lMatrix[computed_i_index][2]
    assert str(ComputedTimePoint) == lMatrix[computed_i_index][3]
    
    # The fold change will be evaluated for intervals [i + x for x in delta_i]
    # delta_i is of length normWindowSize centered on 0
    delta_i = [x for x in range(- normWindowSize // 2 + 1, normWindowSize // 2 + 1)]
    # We need to discard from delta_i all x where i + x is negative or above the len(lRawTimePoint)
    delta_i = [x for x in delta_i if x >=0 and x <len(lRawTimePoint)]
    # Or which corresponds to different chromosome
    delta_i = [x for x in delta_i
               if lMatrix[target_i_index + x][0] == line[0] and
                  lMatrix[computed_i_index + x][0] == line[0]]

    # Sum the FoldChange over each bin:
    sumFoldChange = 0.0
    for x in delta_i:
        # Check they correspond to same interval
        assert lMatrix[target_i_index + x][1] == lMatrix[computed_i_index + x][1]
        assert lMatrix[target_i_index + x][2] == lMatrix[computed_i_index + x][2]

        # If one of the value is too small it is substituted by 1
        if lMatrix[target_i_index + x][4] is None or lMatrix[target_i_index + x][4] < 1:
            lMatrix[target_i_index + x][4] = 1

        if lMatrix[computed_i_index + x][4] is None or lMatrix[computed_i_index + x][4] < 1:
            lMatrix[computed_i_index + x][4] = 1

        sumFoldChange = sumFoldChange + \
            (lMatrix[computed_i_index + x][4] / lMatrix[target_i_index + x][4])

    # Compute the average
    VariationRatio = sumFoldChange / len(delta_i)

    # As we will do a division
    # We cannot get 0
    if VariationRatio == 0.0:
        VariationRatio = 1

    if line[3] == None:
        line[3] = 0

    ModifiedIntervals.append(
        [line[0], line[1], line[2], (line[3] / VariationRatio)])

output_file_name = sFPath_Bigwig[0:len(sFPath_Bigwig) - 7] + \
                       '_HCN_' + str(TargetTimePoint) + 'h.bedgraph'

with open(output_file_name, 'w') as fResult:
    for Interval in ModifiedIntervals:
        fResult.write("%s\n" % (Interval[0] + "\t" +
                                Interval[1] + "\t" + Interval[2] + "\t" + str(Interval[3])))
