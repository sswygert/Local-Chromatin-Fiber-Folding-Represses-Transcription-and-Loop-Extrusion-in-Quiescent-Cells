#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2019 dejun <dejun.lin@gmail.com>
# Created at 2019-09-19 19:25 on dejun@dejun-GS60-2QE
# Usage: contactDecay.py
# Description: Compute contact decay by genes
#
# Distributed under terms of the MIT license.
import os
import numpy as np
from scipy import stats
import math
import sys
import argparse
import subprocess
import re
import pandas as pd
import cooler
import h5py
from pybedtools import BedTool
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10})
from utils import contactDecay


parser = argparse.ArgumentParser()
parser.add_argument("fmcool", type=str,
                    help="Input cooler multiple-binsize contact files")
parser.add_argument("binSize", type=int,
                    help="Use data corresponding to this bin size in fmcool\
                    one of 200, 500, 1000, 3200 and 5000 bp")
parser.add_argument("foutprefix", type=str,
                    help="Output results to files prefixed by this")
parser.add_argument("--bNormNpairs", action='store_true', default=False,
                    help="Whether to normalize the number of contacts by the number of pairs\
                    of each contact matrix diagonal")
parser.add_argument("--fNucSgr", type=str, default="",
                    help="Compute inter-nucleosome contact decay with the nucleosome\
                    peaks defined in this file")
parser.add_argument("--fbed", type=str, default="",
                    help="Optionally compute the contact histogram for genes in this\
                    annotation")
parser.add_argument("--col", type=int, default=-1,
                    help="When --fbed is specified, group contact histogram by this column in the fbed")

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
foutprefix = args.foutprefix
bNormNpairs = args.bNormNpairs
fNucSgr = args.fNucSgr
fbed = args.fbed
col = args.col


assert not (fbed != "" and fNucSgr != ""),\
    "bed-file mode and nucleasome peak mode can't be used at the same time"

# Read the mcool file
mcool = h5py.File(fmcool, 'r')
cool = cooler.Cooler(mcool['resolutions'][str(binSize)])

# save the whole-genome contact decay
nBPMaxWhole = cool.chromsizes.max()
nBinsMaxWhole = math.ceil(nBPMaxWhole / binSize) + 1

headerGroup = header
headerMean = header
headerMedian = header

nucPeaks = None
if fNucSgr != "":
    nucPeaks = pd.read_csv(fNucSgr, sep="\s+", header=None, usecols=(0,1), names=["chr", "loc"])
    # sgr file has 1-based genomic positions which are converted to 0-based for binning 
    # to Cooler bins, which are 0-based and lower bound inclusive and upper bound exclusive
    nucPeaks['loc'] -= 1
    nBinsMaxWhole = nucPeaks.groupby('chr').count().max()

contactDecayWhole = np.zeros(nBinsMaxWhole, dtype=float)

for chromname in cool.chromnames:
    cd = contactDecay(cool, chromname, bNormNpairs=bNormNpairs, nucPeaks=nucPeaks)
    contactDecayWhole[:cd.size] += cd
# trim all the trailing zero elements
contactDecayWhole = np.trim_zeros(contactDecayWhole, trim='b')
contactDecayWholeNormed = contactDecayWhole / contactDecayWhole.sum()
distWhole = np.arange(contactDecayWhole.size)
if fNucSgr == "":
    # when not using nucleosome peaks, set the genomic distance axis in the unit of basepair
    # otherwise leave it to be nucleasome peak index
    distWhole *= binSize
    header += "# GenomicDistance Counts Prob"
else:
    header += "# nucleasomeId Counts Prob"
np.savetxt(foutprefix + ".WholeGenome.contactDecay.mean",
           np.vstack((distWhole, contactDecayWhole, contactDecayWholeNormed)).T,
           "%30.15e", header=header)

if fbed != "":
    # Read the bed file
    bed = BedTool(fbed).to_dataframe()
    # remove '#' from the 1st column name
    bed.columns = bed.columns.str.replace('.*chrom$', 'chrom', regex=True)
    bed['length'] = bed.apply(lambda gene: gene.end - gene.start, axis=1)
    nBPMax = bed['length'].max()
    # plus one because the 1st bin is intra-bin histogram count 
    nBinsMax = math.ceil(nBPMax / binSize) + 1
    # parse number of categories in 'col'
    colName = bed.columns[col]
    colGroups = bed.groupby(colName).size()
    nGroups = colGroups.size
    colGroupsDF = pd.DataFrame({'group': colGroups.index,
                                'nGenes': colGroups.values})
    nHists = colGroups.size
    # save the per-gene contact decay 
    contactDecays = np.zeros((bed.shape[0], nBinsMax), dtype=float)

    # go through each gene and do histogramming
    for index, gene in bed.iterrows():
        key = gene['chrom'] + ':' + str(gene['start']) + '-' + str(gene['end'])
        cd = contactDecay(cool, key, bNormNpairs=bNormNpairs)
        contactDecays[index, :cd.size] = cd
    # normalize the per-gene contact decay
    contactDecayNormed = contactDecays / contactDecays.sum(axis=1)[:, np.newaxis]
    # replace all the empty contact decay with zero
    np.nan_to_num(contactDecayNormed, copy=False)

    distGroup = distWhole[:contactDecays.shape[1]]
    headerGroup += "\n#Row: GenomicDistance: " + str(distGroup).replace('\n', '')
    headerGroup += "\n#Col: genes\n"
    np.savetxt(foutprefix+".contactDecay", contactDecays, "%5d",
               header=headerGroup)
    np.savetxt(foutprefix+".contactDecay.normed", contactDecayNormed, "%30.15e",
               header=headerGroup)

    # plot mean and std
    plt.clf()
    plt.cla()
    # exclude the first bin, i.e., the intra-bin contacts in the plot
    # and limit the number of bins to 100
    nBinsMaxPlot = min(100, nBinsMax - 1)
    binsPlot = np.arange(nBinsMaxPlot)
    headerMean += "# GenomicDistance StdProbOfGroup MeanProbOfGroup\n"
    for i in range(nGroups):
        subsetIndices = bed[colName] == colGroups.index[i]
        contactDecay = contactDecayNormed[subsetIndices].copy()
        meanContactDecay = np.mean(contactDecay, axis=0)
        stdContactDecay = np.std(contactDecay, axis=0)
        eb = plt.errorbar(binsPlot, meanContactDecay[1:(nBinsMaxPlot + 1)],
                          yerr=stdContactDecay[1:(nBinsMaxPlot + 1)],
                          errorevery=5,
                          label=colGroups.index[i])
        eb[-1][0].set_linestyle('--')
        np.savetxt(foutprefix + ".contactDecay." + str(colGroups.index[i]) + ".mean",
                   np.vstack((distGroup, stdContactDecay, meanContactDecay)).T,
                   "%30.15e", header=headerMean)
    plt.xlabel("bins (binSize = %d bp)" % binSize)
    plt.ylabel("Mean contact Probability")
    plt.legend()
    plt.savefig(foutprefix+".contactDecay.ByGroup.Mean.png", dpi=300)

    # plot median and mad
    plt.clf()
    plt.cla()
    headerMedian += "# GenomicDistance MadProbOfGroup MedianProbOfGroup\n"
    for i in range(nGroups):
        subsetIndices = bed[colName] == colGroups.index[i]
        contactDecay = contactDecays[subsetIndices].copy()
        medianContactDecay = np.median(contactDecay, axis=0)
        madContactDecay = stats.median_absolute_deviation(contactDecay, axis=0)
        medianContactDecayNorm = medianContactDecay.sum()
        medianContactDecayPlot = medianContactDecay[1:(nBinsMaxPlot + 1)] / medianContactDecayNorm
        madContactDecayPlot = madContactDecay[1:(nBinsMaxPlot + 1)] / medianContactDecayNorm
        eb = plt.errorbar(binsPlot, medianContactDecayPlot,
                          yerr=madContactDecayPlot,
                          errorevery=5,
                          label=colGroups.index[i])
        eb[-1][0].set_linestyle('--')
        np.savetxt(foutprefix + ".contactDecay." + str(colGroups.index[i]) + ".median",
                   np.vstack((distGroup, madContactDecay, medianContactDecay)).T,
                   "%30.15e", header=headerMedian)
    plt.xlabel("bins (binSize = %d bp)" % binSize)
    plt.ylabel("Median contact probability")
    plt.legend()
    plt.savefig(foutprefix+".contactDecay.ByGroup.Median.png", dpi=300)

    plt.clf()
    plt.cla()
    plt.xlabel("bins (binSize = %d bp)" % binSize)
    plt.ylabel("Genes")
    plt.imshow(contactDecayNormed[:, 1:(nBinsMaxPlot + 1)],
               interpolation='nearest', aspect='auto',
               cmap='plasma', norm=LogNorm(1e-3, 1))
    cb = plt.colorbar()
    cb.set_label("Contact Probability")
    plt.savefig(foutprefix+".contactDecay.normed.png", dpi=300)
