"""
@author: hocinerekaik

Generate coverage matrix of given intervals and time-points. If needed, can complete coverage in missing time-points.

Modified May 5 2022 by lldelisle
To propose command line and dynamic time setting
"""

import numpy as np
import pyBigWig
import re
import tkinter as tk
from tkinter import filedialog, simpledialog
import argparse
import os

argp = argparse.ArgumentParser(
    description=("Generate coverage matrix of given intervals and time-points"))
argp.add_argument('--bigwig', default=None, nargs='+',
                  help="bigwig files (separated by space)")
argp.add_argument('--timepoints', default=None, nargs='+', type=int,
                  help="time points of bigwig files (separated by space)")
argp.add_argument('--bed', default=None,
                  help="bed file with regions on which coverage is computed")
argp.add_argument('--newTimeInterval', default=None, type=int,
                  help="new time interval to which coverage will be computed")
argp.add_argument('--prefix', default="",
                  help="prefix to use for matrix and pdf to build the film")
argp.add_argument('--ini', default=None,
                  help="ini file to use for pyGenomeTracks")
args = argp.parse_args()

# Initialize root
root = None
lFilesList = []

if args.bigwig is None:
    # If no argument is given use interactive mode
    """
    Select bigwig files
    """
    root = tk.Tk()
    root.withdraw()

    lTimePoints = []

    # Let the user enter the time point
    currentTimePoint = simpledialog.askinteger("Time point selection", 'Enter the first time-point in hours')

    while currentTimePoint is not None:
        sFPath_Bigwig = filedialog.askopenfilename(
            title="Select BigWig file " + str(currentTimePoint) + "h")
        root.update()
        if sFPath_Bigwig == "":
            break
        lFilesList.append([sFPath_Bigwig, currentTimePoint, 0.0])
        lTimePoints.append(currentTimePoint)
        currentTimePoint = simpledialog.askinteger("Time point selection",
                                                   'Enter next time-point in hours\nLast one selected = ' +
                                                   str(currentTimePoint) + 'h\nUse cancel when finished')
else:
    lTimePoints = args.timepoints
    if lTimePoints is None:
        raise Exception("--bigwig was set, --timepoints must be set")
    if len(args.bigwig) != len(lTimePoints):
        raise Exception("The number of bigwig file must correspond to the number of timepoints")
    # Check files exists:
    for tp, bw in zip(lTimePoints, args.bigwig):
        if not os.path.exists(bw):
            raise Exception(bw + " does not exist")
        lFilesList.append([bw, tp, 0.0])

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

# Sort the intervals and bigwigs:
lTimePoints, lFilesList = zip(*sorted(zip(lTimePoints, lFilesList)))

"""
compute coverage from Bigwig files
"""
# bed file intervals covering regions of interest
with open(sFPath_Bed) as f:
    lBinsLines = f.read().splitlines()

mTimePoint = []
# mTimePoint is an array
# Each element corresponds to a time-point
# For each time-point, there is an array
# Each element of this array is:
# [chr, start, end, time_point, average]
for sample in lFilesList:
    # sample is [bw_file, time_point, 0.0]
    bw = pyBigWig.open(sample[0])
    result = []
    for line in lBinsLines:
        lSplittedLine = line.split("\t")
        result.append([lSplittedLine[0], lSplittedLine[1], lSplittedLine[2], sample[1], bw.stats(
            lSplittedLine[0], int(lSplittedLine[1]), int(lSplittedLine[2]), type="mean")[0]])
    mTimePoint.append(result)
    bw.close()


"""
complete missing Time-points
"""
if args.newTimeInterval is None:
    if root is None:
        root = tk.Tk()
        root.withdraw()
    newTimeInterval = simpledialog.askinteger("Time interval selection", 'Enter the desired new time interval in hours')
    root.update()
else:
    newTimeInterval = args.newTimeInterval

# Store all time points computed:
mTimePointComputed = []
# Store the last result
# in TimePoint
previousMasterTimepointResults = None
previousMasterTimepoint = None
for i, TimePoint in enumerate(mTimePoint):
    # TimePoint is an array corresponding to a single time_point,
    # each element is like [chr, start, end, time_point, average]
    current_timepoint = TimePoint[0][3]
    # For the first one we simply store the TimePoint
    if current_timepoint == lTimePoints[0]:
        mTimePointComputed.append(TimePoint)
        previousMasterTimepointResults = TimePoint
        previousMasterTimepoint = current_timepoint
        continue
    # For the other we will compute the missing timepoints
    # Between the previousMasterTimepoint + 1 and
    # current_timepoint (included)
    for nTimingIndex, vActualTime in enumerate(np.arange(previousMasterTimepoint + newTimeInterval,
                                                         current_timepoint + 1,
                                                         newTimeInterval),
                                               start=1):
        lNewTimePoint = []
        # For each interval
        for nInterval, line in enumerate(TimePoint):
            # If the value for the current_timepoint is None we put 0
            if line[4] == None:
                line[4] = 0
            # If the value for the previousMasterTimepointResult is None we put 0
            if previousMasterTimepointResults[nInterval][4] == None:
                previousMasterTimepointResults[nInterval][4] = 0.0
            # We get the delta coverage per one newTimeInterval
            # (Previous value - current value) / ((time current - time previous) / newTimeInterval)
            vDeltaCoverage = (previousMasterTimepointResults[nInterval][4] - line[4]) / ((current_timepoint - previousMasterTimepoint) / newTimeInterval)
            # We append the new value:
            # previous value - (delta * nTimingIndex)
            lNewTimePoint.append(
                [line[0], line[1], line[2], vActualTime, previousMasterTimepointResults[nInterval][4] - (vDeltaCoverage * nTimingIndex)])
        # Once new values are computed for all intervals we add it to mTimePointComputed
        mTimePointComputed.append(lNewTimePoint)
    # Update the previous master:
    previousMasterTimepointResults = TimePoint
    previousMasterTimepoint = current_timepoint


"""
Generate Bedgraph files for pyGenomeTracks figures and save matrix
"""

dirname = os.path.dirname(args.prefix)
if dirname != '' and not os.path.isdir(dirname):
    print(dirname + " does not exists. Will write to current directory.")
    args.prefix = os.path.basename(args.prefix)

fMatrix = open(args.prefix + 'HoxAll_' +
               str(newTimeInterval) + 'h_timepoint.txt', 'w')
for TimePoint in mTimePointComputed:
    # TimePoint is an array corresponding to a single time_point,
    # each element is like [chr, start, end, time_point, average]
    if args.ini is not None:
        # Write the bedgraph to current.bedgraph
        fResult = open('current.bedgraph', 'w')
    for Interval in TimePoint:
        if args.ini is not None:        
            fResult.write("%s\n" % (Interval[0] + "\t" +
                                    Interval[1] + "\t" + Interval[2] + "\t" + str(Interval[4])))
        fMatrix.write("%s\n" % (Interval[0] + "\t" +
                                Interval[1] + "\t" + Interval[2] + "\t" + str(Interval[3]) + "\t" + str(Interval[4])))
    if args.ini is not None:
        fResult.close()
        os.system("pyGenomeTracks --tracks " + args.ini + " --region chr2:74,664,307-74,771,562 --outFileName " +
                  args.prefix + "_" + str(Interval[3]) + ".png --title " + str(Interval[3]) + "h")
fMatrix.close()
