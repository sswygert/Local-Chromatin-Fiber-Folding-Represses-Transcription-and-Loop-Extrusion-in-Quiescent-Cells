#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2020 dejun <dejun.lin@gmail.com>
# Created at 2020-09-14 17:51 on dejun@dejun-GS60-2QE
# Usage: pileupContacts.py
# Description: Pileup contact from Cooler contact matrix of genomic loci
# specified in a bed file
#
# Distributed under terms of the MIT license.
import os
import scipy as sp
import numpy as np
import math
import sys
import argparse
import warnings
import subprocess
import re
from pybedtools import BedTool
import pandas as pd
import cooler
import h5py
import matplotlib
import matplotlib.colors as colors
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from utils import readMcool, pixels2Coo, cool2pixels, getSubCoo
from scipy.sparse import coo_matrix, csr_matrix, csc_matrix

parser = argparse.ArgumentParser()
parser.add_argument("fcool", type=str,
                    help="Input Cooler file of the contact matrix")
parser.add_argument("fbed", type=str,
                    help="Input bed file of a list of genomic loci\
                    whose contacts are to be piled up for output")
parser.add_argument("nBins", type=int,
                    help="The output will be a contact matrix of nBins x nBins\
                    centering on each peak from the bed file")
parser.add_argument("foutput", type=str,
                    help="Output the piled-up contact matrix to this file")
parser.add_argument("--binSize", type=int, default=-1,
                    help="Use this to select the bin size from the input mcool\
                    file. Default to -1, meaning that the inputs are treated as\
                    single-binsize .cool files")
parser.add_argument("--title", type=str, default=None,
                    help="Set title of the plot")
parser.add_argument("--bMedian", action='store_true', default=False,
                    help="If set, output the median of the submatrices;\
                    otherwise output the mean")
parser.add_argument("--gDistNBins", type=int, default=0,
                    help="Only include the submatrices whose center is equal or below\
                    this number of bins from the major diagonal")
parser.add_argument("--excludeZeroSub", type=float, default=0,
                    help="If set, exclude submatrices from those that pass the\
                    --gDistNBins threshold which have less than this\
                    percent zero entries of the nBins*nBins total number of\
                    entries. A value of zero means all submatrices are included")
parser.add_argument("--bNormDiag", action='store_true', default=False,
                    help="If set, the submatrices are normalized by subtracting the\
                    median and then dividing by the median absolute deviation\
                    of their  respective diagonal. This has the effect\
                    of subtracting the background close to the diagonal")
parser.add_argument("--bIntra", action='store_true', default=False,
                    help="If set, the submatrices corresponding to the intra\
                    loci will be included in the pileup")
parser.add_argument("--bInter", action='store_true', default=False,
                    help="If set, the submatrices corresponding to the inter\
                    loci will be included in the pileup")
parser.add_argument("--cBounds", nargs=2, type=float, default=None,
                    help="Min and max for the output color bar")
parser.add_argument("--cLog", action='store_true', default=False,
                    help="If set, use log scale in the output color bar")

args = parser.parse_args()

header = "#"+" ".join(sys.argv)+"\n"

# Check if current script is under revision control
gitls = subprocess.Popen('cd '+os.path.dirname(os.path.realpath(__file__)) +
                         ' && git ls-files --error-unmatch ' + os.path.realpath(__file__),
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
fbed = args.fbed
nBins = args.nBins
foutput = args.foutput
binSize = args.binSize
title = args.title
bMedian = args.bMedian
gDistNBins = args.gDistNBins
excludeZeroSub = args.excludeZeroSub
bNormDiag = args.bNormDiag
bIntra = args.bIntra
bInter = args.bInter
cBounds = args.cBounds
cLog = args.cLog

assert bIntra or bInter, "Both intra- and inter-loci pileups are off"

cool, binSize = readMcool(fcool, binSize)
bins = cool.bins()
matrices = cool.matrix(balance=False, sparse=True)
p = cool2pixels(cool)
bed = BedTool(fbed).to_dataframe()
# set up the piled-up contact matrix for output
nBps = nBins * binSize
nBinsZeroExclude = int(nBins * nBins * excludeZeroSub)
# report how many submatrices actually contribute to the
# pileup matrix
nSubs = 0
# loop over each chromosome and perform intra-chromosomal pileup
msOut = []
nChrs = len(cool.chromnames)
for iChr in range(nChrs):
    iChrName = cool.chromnames[iChr]
    iBed = bed[bed.chrom == iChrName]
    # the double loop below requires the entries in the bed file
    # to be sorted because Cooler file usually stores the upper triangle
    # of the contact matrix. For inter chromosomal pileup, we also need to
    # make sure the chromosome are sorted too
    assert np.all(np.diff(iBed['start']) >= 0),\
        f"Entries in {fbed} is not sorted"
    if bNormDiag:
        # precompute the background contacts at each diagonal
        iPixels = p.fetch(iChrName)
        iNBinsTotal = bins.fetch(iChrName).shape[0]
        iPixels['binDist'] = iPixels['bin2_id'] - iPixels['bin1_id']
        diagGroups = iPixels.groupby('binDist')['count']
        # to avoid missing diagonal sums, need to initialize the whole array
        normDiags = np.zeros(iNBinsTotal, dtype=float)
        # use median or mean as the background depending on which to output
        if bMedian:
            dfNormDiags = diagGroups.median()
        else:
            dfNormDiags = diagGroups.mean()
        normDiags[dfNormDiags.index] = dfNormDiags
    # build the region string first
    regions = []
    for index, gene in iBed.iterrows():
        iCenter = (gene['end'] + gene['start']) // 2
        # deal with the edges
        if iCenter - nBps // 2 <= 0:
            iStart = 0
            iLast = iStart + nBps
        elif iCenter + nBps // 2 > cool.chromsizes[iChrName]:
            iLast = cool.chromsizes[iChrName]
            iStart = iLast - nBps + (iLast % binSize != 0) * binSize
        else:
            iStart = iCenter - nBps // 2
            iLast = iStart + nBps - (iStart % binSize != 0) * binSize
        region = gene['chrom'] + ':' + str(int(iStart)) + '-' + str(int(iLast))
        regions.append(region)
        assert bins.fetch(regions[-1]).shape[0] == nBins,\
            f"Sub contact matrix for {gene.values} in {fbed} has wrong shape:"\
            f" expecting {nBins} but got {bins.fetch(regions[-1]).shape[0]}"
    # loop over pairs of submatrices and aggregate
    for i in range(len(regions) - 1):
        iRegion = regions[i]
        iBins = bins.fetch(iRegion).index.values
        iBinCenter = iBins[nBins // 2]
        jMin = i if bIntra else i + 1
        jMax = len(regions) if bInter else i + 1
        for j in range(jMin, jMax):
            jRegion = regions[j]
            jBins = bins.fetch(jRegion).index.values
            jBinCenter = jBins[nBins // 2]
            binDist = abs(jBinCenter - iBinCenter)
            if binDist > gDistNBins:
                continue
            mSub = matrices.fetch(iRegion, jRegion)
            # exclude the all zero submatrices
            if mSub.nnz < nBinsZeroExclude:
                continue
            if bNormDiag:
                # normalize each entry in mSub by the respective median of their
                # corresponding diagonal
                iDiags = -np.subtract.outer(iBins, jBins)
                msOut.append(mSub.toarray() - normDiags[iDiags])
            else:
                msOut.append(mSub.toarray())
            nSubs += 1
msOut = np.array(msOut)

if bMedian:
    mOut = np.median(msOut, axis=0)
    cbLabel = "Median "
else:
    mOut = np.mean(msOut, axis=0)
    cbLabel = "Mean "

if bNormDiag:
    cbLabel += "normalized counts"
else:
    cbLabel += "counts"

# save as dense matrix
header += f"# Number of submatrices included = {nSubs}\n"
np.savetxt(foutput, mOut, "%30.15e", header=header)

plt.rcParams.update({'font.size': 22})
fig = plt.figure(figsize=(13,10),dpi=200)
ax = plt.gca()

if cBounds is not None:
    vMin = cBounds[0]
    vMax = cBounds[1]
else:
    vMin = mOut.min()
    vMax = mOut.max()

if cLog:
    norm = colors.LogNorm(vmin=vMin, vmax=vMax)
else:
    norm = None

plt.rcParams.update({'font.size': 22})
fig = plt.figure(figsize=(13,10),dpi=200)
ax = plt.gca()

sns.heatmap(mOut,
            vmin=vMin, vmax=vMax, norm=norm,
            cmap='plasma', cbar_kws={'label' : cbLabel},
            ax=ax,
            linewidths=0.0, rasterized=True
            )

# use kb as unit of ticks and only plot a tick evyer 2 kb
every = 2000 // binSize
tickLabels = ((np.arange(nBins) - nBins // 2) * binSize / 1000)[every::every]
tickPositions = np.arange(nBins)[every::every]
ax.set_xticks(ticks=tickPositions)
ax.set_xticklabels(labels=tickLabels)
ax.set_yticks(ticks=tickPositions)
ax.set_yticklabels(labels=tickLabels)

ax.set_xlabel(f"Genomic distance (kb)")
ax.set_ylabel(f"Genomic distance (kb)")
ax.set_title(title)
fig.savefig(foutput + ".svg", dpi=300, metadata={"Title" : header})
