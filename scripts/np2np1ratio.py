#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2020 dejun <dejun.lin@gmail.com>
# Created at 2020-03-23 11:05 on dejun@dejun-GS60-2QE
# Usage: np2np1ratio.py
# Description: ratio of number of contacts at N+2 and N+1
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
import pandas as pd
import cooler
import h5py
import matplotlib
matplotlib.use('Agg')
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10})
from scipy.sparse import coo_matrix, csr_matrix, lil_matrix
import pyBigWig
from utils import readMcool, cool2pixels, getSubCoo

parser = argparse.ArgumentParser()
parser.add_argument("fmcool", type=str,
                    help="Input cooler multiple-binsize contact files")
parser.add_argument("binSize", type=int,
                    help="Use data corresponding to this bin size in fmcool\
                    one of 10, 100, 150, 200, 500, 1000, 3200 and 5000 bp")
parser.add_argument("windowSize", type=int,
                    help="Bin the counts by window of this number of bp")
parser.add_argument("foutprefix", type=str,
                    help="Output results to files prefixed by this")
parser.add_argument("--fNucSgr", type=str, default="",
                    help="Sgr file for the location of the nucleosome")
parser.add_argument("--np2np1Bins", type=int, nargs=4, default=[0,0,0,0],
                    help="Genomic distance ranges for n+1 and n+2 nucleosome\
                    contacts. The arguments are:\
                    [np1bpBegin, np1bpLast,  np2bpBegin, np2bpLast]")

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

fmcool = args.fmcool
binSize = args.binSize
windowSize = args.windowSize
foutprefix = args.foutprefix
fNucSgr = args.fNucSgr
np2np1Bins = args.np2np1Bins

assert not (fNucSgr == "" and np2np1Bins == [0,0,0,0]),\
    "Provide either a Sgr nucleosome peak file or the definition\
    of nucleosome bins"

assert not (fNucSgr != "" and np2np1Bins != [0,0,0,0]),\
    "Can't work with both a Sgr nucleosome peak file\
    and definition of nucleosome bins. Use only one\
    of these two."

assert binSize <= 200, "binSize is larger than 200 bp"
assert windowSize >= binSize, "windowSize is smaller than binSize"

nBinsPerWindow = windowSize // binSize

cool, binSize = readMcool(fmcool, binSize)
p = cool2pixels(cool)
bins = cool.bins()
nBins = bins[:].groupby('chrom').size()
# store the results into a sparse matrix
chromnames = list(filter(re.compile(r"[^M]").search, cool.chromnames))
# rows are chromosome; cols are genomic position in bp
# I don't use bin index because outputing to bigwig needs genomic
# position
row1, pos1, data1 = [], [], []
row2, pos2, data2 = [], [], []

nucPeaks = None
if fNucSgr != "":
    nucPeaks = pd.read_csv(fNucSgr, sep="\s+", header=None, usecols=(0,1), names=["chr", "loc"])
    # assign nucleosome peaks to bins
    # allow +/- 100 bp around the peak to be called a nucleosome
    nucPeaks['binMin'] = (nucPeaks['loc'] - 100) // binSize
    nucPeaks['binMax'] = (nucPeaks['loc'] + 100) // binSize
    nucPeaks.loc[nucPeaks['binMin'] < 0, 'binMin'] = 0
    # prevent neighboring bins from overlapping
    def trimBinEdges(df : pd.DataFrame):
        df.loc[df['binMax'].shift(1) >= df['binMin'], 'binMin'] =\
            df['binMax'].shift(1, fill_value=0) + 1
        return df
    nucPeaks = nucPeaks.groupby('chr').apply(trimBinEdges)
elif np2np1Bins != [0,0,0,0]:
    np1bpBegin, np1bpLast,  np2bpBegin, np2bpLast = np2np1Bins


for iChr in range(len(chromnames)):

    chromname = chromnames[iChr]
    sp = p.fetch(chromname)
    sb = bins.fetch(chromname)
    nWindows = sb.shape[0] - nBinsPerWindow + 1
    sp['binDist'] = sp.apply(lambda c: c.bin2_id - c.bin1_id, axis=1)
    # index by genomic distance
    spi = sp.set_index('binDist')

    if nucPeaks is not None:
        # assign the loci pair to nucleosome index
        iNuc = nucPeaks.loc[nucPeaks['chr'] == "chr"+chromname, ['binMin', 'binMax']]. \
            reset_index(drop=True)
        # offset the nucleosome bin index to be consistent with the cooler file
        iNuc += sb.index[0]
        nucBins = pd.IntervalIndex.from_arrays(iNuc['binMin'], iNuc['binMax'], closed='both')
        sp['bin1NucId'] = pd.cut(sp['bin1_id'], nucBins, labels=False).\
            apply(lambda df : nucBins.get_loc(df))
        sp['bin2NucId'] = pd.cut(sp['bin2_id'], nucBins, labels=False).\
            apply(lambda df : nucBins.get_loc(df))
        # drop loci pairs that aren't assigned to a nucleosome
        sp = sp.loc[~pd.isna(sp['bin1NucId']) & ~pd.isna(sp['bin2NucId'])].reset_index(drop=True)
        sp['bin1NucId'] = sp['bin1NucId'].astype(int)
        sp['bin2NucId'] = sp['bin2NucId'].astype(int)
        # assign the contacts to either n+1 or n+2
        # and sum the total number of n+1 and n+2 contacts for each loci
        np1 = sp.loc[sp['bin2NucId'] == sp['bin1NucId'] + 1, ['bin1_id', 'count']].\
            reset_index(drop=True).\
            groupby('bin1_id').sum()
        np2 = sp.loc[sp['bin2NucId'] == sp['bin1NucId'] + 2, ['bin1_id', 'count']].\
            reset_index(drop=True).\
            groupby('bin1_id').sum()
    else:
        # select the contacts that are [np1bpBegin, np1bpLast] bp away
        np1 = spi.loc[spi.index.isin(range(max(1, np1bpBegin//binSize), np1bpLast//binSize+1))]. \
          reset_index(). \
          drop(columns=['bin2_id', 'binDist']). \
          groupby('bin1_id').sum()
        # select the contacts that are [np2bpBegin, np2bpLast] bp away
        np2 = spi.loc[spi.index.isin(range(max(2, np2bpBegin//binSize), np2bpLast//binSize+1))]. \
          reset_index(). \
          drop(columns=['bin2_id', 'binDist']). \
          groupby('bin1_id').sum()

    # keep the starting genomic position of each bin
    np1['start'] = sb.loc[np1.index, 'start']
    np2['start'] = sb.loc[np2.index, 'start']
    # reset bin index so that it always starts with zero for each chrom
    np1 = np1.set_index(np1.index - sb.index[0])
    np2 = np2.set_index(np2.index - sb.index[0])
    # histogram the contacts into the windows
    for iWindowOffset in range(nBinsPerWindow):
        # each loop step contributes to iBin - iWindowOffset windows
        iNp1 = np1.loc[(np1.index - iWindowOffset).isin(range(0, nWindows))]
        iCol1 = iNp1.index
        iRow1 = np.repeat(iChr, iCol1.size)
        iData1 = iNp1['count'].to_numpy()
        iPos1 = iNp1['start'].to_numpy()
        row1.extend(iRow1)
        data1.extend(iData1)
        pos1.extend(iPos1)
        iNp2 = np2.loc[(np2.index - iWindowOffset).isin(range(0, nWindows))]
        iCol2 = iNp2.index
        iRow2 = np.repeat(iChr, iCol2.size)
        iData2 = iNp2['count'].to_numpy()
        iPos2 = iNp2['start'].to_numpy()
        row2.extend(iRow2)
        data2.extend(iData2)
        pos2.extend(iPos2)

p1 = csr_matrix((data1, (row1, pos1)), shape=(len(chromnames), cool.chromsizes.max() - windowSize + 1))
p2 = csr_matrix((data2, (row2, pos2)), shape=(len(chromnames), cool.chromsizes.max() - windowSize + 1))
ratio = np.array(np.nan_to_num(p2 / p1, nan=0.0, posinf=0.0, neginf=0.0))
ratioMean = ratio / ratio[np.nonzero(ratio)].mean()
p12 = np.array(p2 - p1)
p1p2 = np.array(p2 + p1)
ratioDiffSum = np.array(np.nan_to_num(p12 / p1p2, nan=0.0, posinf=0.0, neginf=0.0))

# output to bigwig file
p1Coo = p1.tocoo()
p2Coo = p2.tocoo()
ratioCoo = coo_matrix(ratio)
ratioMeanCoo = coo_matrix(ratioMean)
ratioDiffSumCoo = coo_matrix(ratioDiffSum)

def coo2bigwig(mCoo : coo_matrix, chrNameSizes : pd.DataFrame, binSize : int, fname : str):
    """Create a bigwig file from a coo matrix

    Args:
        mCoo: `coo_matrix` Input coo_matrix with each row corresponding
        to one chromosome and each column to one genomic position
        chrNameSizes: `pd.DataFrame` with two columns: name and size of
        the chromosomes
        binSize: binSize of the contact matrix
        fname: output bigwig to this file

    Returns:
    """
    chrs = chrNameSizes['name'].to_numpy().astype(str)[mCoo.row]
    iSorted = np.argsort(chrs, kind='mergesort')
    chrs = chrs[iSorted]
    starts = mCoo.col[iSorted].astype(np.int64)
    # can't use overlapping windows for bigwig file
    # so pretend the window to be non-overlapping
    ends=(starts + binSize).astype(np.int64)
    values=mCoo.data[iSorted].astype(float)
    bw = pyBigWig.open(fname, "w")
    bw.addHeader(list(chrNameSizes.itertuples(index=False, name=None)))
    bw.addEntries(chrs, starts, ends=ends, values=values)
    bw.close()

coo2bigwig(p1Coo, cool.chroms()[:], binSize, foutprefix+".np1.bw")
coo2bigwig(p2Coo, cool.chroms()[:], binSize, foutprefix+".np2.bw")
coo2bigwig(ratioCoo, cool.chroms()[:], binSize, foutprefix+".np2np1.bw")
coo2bigwig(ratioMeanCoo, cool.chroms()[:], binSize, foutprefix+".np2np1normed.bw")
coo2bigwig(ratioDiffSumCoo, cool.chroms()[:], binSize, foutprefix+".np2np1diffsum.bw")
