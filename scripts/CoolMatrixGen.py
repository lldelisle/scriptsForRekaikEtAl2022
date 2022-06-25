"""
@author: hocinerekaik

Generates a balanced cool matrix of a given interval and for each time-point.
Data of missing time-points are generated and plotted with pyGenomeTracks.

Modified May 5 2022 by lldelisle
To propose command line and dynamic time setting
"""

import cooler
from numpy.core.numeric import NaN
import numpy as np
import pandas as pd
import argparse
import os

argp = argparse.ArgumentParser(
    description=("Generate interpolated cool matrices from existing time-points"))
argp.add_argument('--cool', default=None, nargs='+',
                  required=True,
                  help="cool files (separated by space)")
argp.add_argument('--timepoints', default=None, nargs='+', type=int,
                  required=True,
                  help="time points of cool files (separated by space)")
argp.add_argument('--newTimeInterval', default=None, type=int,
                  required=True,
                  help="new time interval to which coverage will be computed")
argp.add_argument('--prefix', default="",
                  help="prefix to use for matrix and pdf to build the film")
argp.add_argument('--region', default=None, required=True,
                  help="region on which to restrict computation")
argp.add_argument('--ini', default=None, required=True,
                  help="ini file to use for pyGenomeTracks")
argp.add_argument('--plottedRegion', default=None, required=True,
                  help="region to use for pyGenomeTracks")
argp.add_argument('--keepCopy', action="store_true",
                  help="keep copy of each cool file")
args = argp.parse_args()

# Check the given args:
if len(args.cool) < 2:
    raise Exception("Requires at least 2 cool files")
if len(args.cool) != len(args.timepoints):
    raise Exception("The number of bigwig file must correspond to the number of timepoints")
# Check files exists:
for cool_fn in args.cool:
    if not os.path.exists(cool_fn):
        raise Exception(cool_fn + " does not exist")

# Sort the intervals and cool files:
lTimePoints, lFilesList = zip(*sorted(zip(args.timepoints, args.cool)))

for i, time_start in enumerate(lTimePoints[:-1]):
    time_end = lTimePoints[i + 1]
    """
    Open balanced cool files
    """
    a = cooler.Cooler(lFilesList[i])
    b = cooler.Cooler(lFilesList[i + 1])

    """""
    Get balanced values
    """""

    a_Values = a.matrix(as_pixels=True, join=True).fetch(args.region)
    b_Values = b.matrix(as_pixels=True, join=True).fetch(args.region)
    ab_Values = pd.merge(a_Values, b_Values, on=['chrom1', 'start1', 'chrom2', 'start2'])
    assert ab_Values['end1_x'].equals(ab_Values['end1_y'])
    assert ab_Values['end2_x'].equals(ab_Values['end2_y'])

    ab_Values['balanced_diff'] = ab_Values['balanced_y'] - ab_Values['balanced_x']

    # Index the bins:
    all_bins_1 = ab_Values[['chrom1', 'start1', 'end1_x']]
    all_bins_2 = ab_Values[['chrom2', 'start2', 'end2_x']]
    all_bins_1.columns = ['chrom', 'start', 'end']
    all_bins_2.columns = ['chrom', 'start', 'end']
    all_bins = pd.concat([all_bins_1, all_bins_2])
    # I add a bin from chrom start as it raises some issues
    all_bins = pd.concat([all_bins,
                          pd.DataFrame({'chrom': [args.region.split(':')[0]],
                                        'start': [0],
                                        'end': min(all_bins['start'])}
                                       )])
    all_bins = all_bins.drop_duplicates().sort_values(['chrom', 'start'])
    all_bins.reset_index(inplace=True,drop=True)

    temp_bins1 = all_bins.copy()
    temp_bins1.columns = ['chrom1', 'start1', 'end1_x']
    temp_bins1['bin1_id'] = temp_bins1.index
    ab_Values = pd.merge(ab_Values, temp_bins1)
    temp_bins1.columns = ['chrom2', 'start2', 'end2_x', 'bin2_id']
    ab_Values = pd.merge(ab_Values, temp_bins1)

    for new_time in np.arange(time_start,
                              time_end,
                              args.newTimeInterval):
        pixels = ab_Values.loc[:, ('bin1_id', 'bin2_id')]
        pixels['count'] = ab_Values['balanced_x'] + ab_Values['balanced_diff'] * (new_time - time_start) / (time_end - time_start)

        cooler.create_cooler('current.cool', all_bins, pixels.dropna(), dtypes = {'count': np.float64})
        os.system("pyGenomeTracks --tracks " + args.ini + " --region " + args.plottedRegion + " --outFileName " +
                  args.prefix + "_" + str(new_time) + ".png --title " + str(new_time) + "h")
        
        if args.keepCopy:
            cooler.create_cooler(args.prefix + "_" + str(new_time) + ".cool", all_bins, pixels.dropna(), dtypes = {'count': np.float64})
