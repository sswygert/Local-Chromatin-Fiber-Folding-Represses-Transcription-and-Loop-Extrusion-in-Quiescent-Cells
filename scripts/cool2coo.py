#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2019 dejunlin <dejun.lin@gmail.com>
# Created at 2019-09-13 14:09 on dejunlin@threonine.gs.washington.edu
# Usage: cool2coo.py
# Description: Convert single-binsize cooler file to coo matrix
#
# Distributed under terms of the MIT license.
import os
import numpy as np
import sys
import argparse
import subprocess
import re
import cooler
import csv
from scipy.sparse import coo_matrix
from scipy.stats import median_abs_deviation
from utils import readMcool, cool2pixels, getSubCoo

parser = argparse.ArgumentParser()
parser.add_argument("fcool", type=str,
                    help="Input cooler file")
parser.add_argument("fcoo", type=str,
                    help="Output coo file")
parser.add_argument("--sub", type=str, default="",
                    help="Subset the contact matrix to this genomic region")
parser.add_argument("--filterMaxMad", type=int, default=0,
                    help="If set, filter out the bins that fall below this number\
                    of median absolute deviations from the median log marginal contacts\
                    in the individual chromosomes. The excluded bins will be \
                    reported as a list of bin indices")
parser.add_argument("--bBalance", action='store_true', default=False,
                    help="If apply matrix balancing (or ICE normalization) to output")
parser.add_argument("--bNormByMax", action='store_true', default=False,
                    help="If divide the contact counts by their maximum")

args = parser.parse_args()

header = "#"+" ".join(sys.argv)+"\n"

# Check if current script is under revision control
gitls = subprocess.Popen('cd '+os.path.dirname(os.path.realpath(__file__)) +
                         ' && git ls-files --error-unmatch ' + __file__,
                         shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, encoding='utf-8')
if not gitls.stderr.read():
    gitrev = subprocess.Popen('cd '+os.path.dirname(os.path.realpath(__file__))
                              + ' && git rev-parse HEAD --abbrev-ref HEAD',
                              shell=True, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, encoding='utf-8')
    if not gitrev.stderr.read():
        p = re.compile("\n")
        header += "# @rev " + p.sub("\n# @branch ",
                                    gitrev.stdout.read().strip()) + "\n"

fcool = args.fcool
fcoo = args.fcoo
sub = args.sub
bBalance = args.bBalance
bNormByMax = args.bNormByMax
filterMaxMad = args.filterMaxMad

assert not (bBalance and bNormByMax),\
    f"Can't use --bBalance and --bNormByMax together\n"

cool, binSize = readMcool(fcool, -1)
pixels = cool2pixels(cool)

if sub == "":
    df = pixels[:]
    shape = pixels.shape
else:
    df = pixels.fetch(sub)
    binsRangeSub = cool.extent(sub)
    nBinsSub = binsRangeSub[1] - binsRangeSub[0]
    shape = (nBinsSub, nBinsSub)

assert not (bBalance and (filterMaxMad> 0) ),\
    "Can't set bBalance and filterMaxMad at the same time"

if filterMaxMad> 0:
    # This is based on the paper https://www.nature.com/articles/nature24281
    # and the cooler.balance_cooler() in
    # https://github.com/open2c/cooler/blob/d6a7d19e5933d5a2a8a31a93bb13c0f06ccb075d/cooler/balance.py#L389-L399
    # Filtering by a high pass filter for number of median absolute deviation
    # This filtering involves genome wide median of the log marginal counts 
    dfAll = pixels[:]
    # Get the marginal counts
    countsMarginal = \
        np.bincount(dfAll['bin1_id'], weights=dfAll['count'],
                    minlength=shape[0]) +\
        np.bincount(dfAll['bin2_id'], weights=dfAll['count'],
                    minlength=shape[0])
    # Normalize each chromosomes' marginal counts by the respective chromosome's
    # median marginal count
    binRangeChroms = cool._load_dset("indexes/chrom_offset")
    for binStart, binEnd in zip(binRangeChroms[:-1], binRangeChroms[1:]):
        countsMarginalChrom = countsMarginal[binStart:binEnd]
        countsMarginal[binStart:binEnd] /= np.median(
            countsMarginalChrom[countsMarginalChrom > 0])
    # Compute the cutoff based on genome wide median of the log marginal counts
    logCountsMarginal = np.log(countsMarginal[countsMarginal > 0])
    medLogCountsMarginal = np.median(logCountsMarginal)
    madLogCountsMarginal = median_abs_deviation(logCountsMarginal)
    cutoff = np.exp(medLogCountsMarginal - filterMaxMad * madLogCountsMarginal)
    binsExcl = np.where(countsMarginal < cutoff)[0]
    binsExclDict = dict.fromkeys(binsExcl)
    # remove the rows from df
    rowsIncl = np.array([ not (i in binsExclDict) for i in df['bin1_id'] ]) &\
        np.array([ not (j in binsExclDict) for j in df['bin2_id'] ])
    df = df[rowsIncl]
    if sub == "":
        binsExclSub = binsExcl
    else:
        # Need to subtract the bin offset index
        binsExclSub = binsExcl[np.where((binsRangeSub[0] <= binsExcl) &
                                        (binsExcl < binsRangeSub[1]))] -\
            binsRangeSub[0]
    header += f"# binsExcluded = {binsExclSub}".replace("\n", "") + "\n"

if bBalance:
    # the parameters were taken from Seungsoo's pipeline as reported on 
    # https://doi.org/10.1016/j.molcel.2018.11.020
    bias, stats = cooler.balance_cooler(cool, mad_max=9, ignore_diags=1)
    # Cooler fill in NaN for bins that are filtered out according to
    # https://github.com/mirnylab/cooler/blob/843dadca5ef58e3b794dbaf23430082c9a634532/cooler/balance.py#L333
    # To prevent propagating NaN downstream, I replace them with zero here
    np.nan_to_num(bias, copy=False)
    df['count'] = bias[df['bin1_id']] * bias[df['bin2_id']] * df['count']
elif bNormByMax:
    df['count'] /= df['count'].max()

header += f"# shape = {shape}\n"

with open(fcoo, 'w') as fh:
    fh.write(header)
    # set bin offset before printing out
    if sub != "":
        df['bin1_id'] -= binsRangeSub[0]
        df['bin2_id'] -= binsRangeSub[0]
    df.to_csv(fh, float_format="%30.15e", quoting=csv.QUOTE_NONE, header=False,
              index=False, sep=' ', escapechar=' ')

