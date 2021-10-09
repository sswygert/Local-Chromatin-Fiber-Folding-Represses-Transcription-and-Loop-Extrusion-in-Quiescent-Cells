#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2019 dejunlin <dejun.lin@gmail.com>
# Created at 2019-11-12 15:04 on dejunlin@threonine.gs.washington.edu
# Usage: groupContacts.py
# Description: Group the contact matrices for regions in an annotation
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
from pybedtools import BedTool
from scipy.sparse import coo_matrix
import matplotlib
matplotlib.use('Agg')
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 8})
from utils import contactDecay

parser = argparse.ArgumentParser()
parser.add_argument("fmcool", type=str,
                    help="Input cooler multiple-binsize contact files")
parser.add_argument("binSize", type=int,
                    help="Use data corresponding to this bin size in fmcool\
                    one of 10, 50, 100, 200, 500, 1000, 3200 and 5000 bp")
parser.add_argument("fbed", type=str,
                    help="Group the contact matrices by genes in this\
                    annotation")
parser.add_argument("every", type=int,
                    help="Group every this number of genes")
parser.add_argument("foutprefix", type=str,
                    help="Output results to files prefixed by this")

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
fbed = args.fbed
every = args.every
foutprefix = args.foutprefix

# Read the mcool file
mcool = h5py.File(fmcool, 'r')
cool = cooler.Cooler(mcool['resolutions'][str(binSize)])
# Read the bed file
bed = BedTool(fbed).to_dataframe()
# remove '#' from the 1st column name
bed.columns = bed.columns.str.replace('.*chrom$', 'chrom', regex=True)
bed['length'] = bed.apply(lambda gene: gene.end - gene.start, axis=1)

def getContactMatrix(cool, region):
    m = cool.matrix(balance=False, sparse=True).fetch(region)
    return m

nGroups = math.ceil(bed.shape[0] / every)
mGroups = []
mGroupsNormed = []
cdGroups = []
cdGroupsNormed = []
nBinsGroups = np.zeros(nGroups, dtype=int)
for iGroup in range(nGroups):
    subBed = bed[iGroup * every : (iGroup + 1) * every]
    # plus one because the 1st bin is intra-bin histogram count 
    nBPMax = subBed['length'].max()
    nBinsMax = math.ceil(nBPMax / binSize) + 1
    mGroup = coo_matrix((nBinsMax, nBinsMax), dtype=int)
    cdGroup = np.zeros(nBinsMax, dtype=float)
    nBinsGroups[iGroup] = nBinsMax
    for index, gene in subBed.iterrows():
        region = gene['chrom'] + ':' + str(gene['start']) + '-' + str(gene['end'])
        m = getContactMatrix(cool, region)
        mGroup.col = np.hstack((mGroup.col, m.col))
        mGroup.row = np.hstack((mGroup.row, m.row))
        mGroup.data = np.hstack((mGroup.data, m.data))
        cd = contactDecay(cool, region)
        cdGroup[:cd.size] += cd
    mGroup.sum_duplicates()
    mGroupNormed = coo_matrix((mGroup.data / mGroup.data.sum(), (mGroup.row, mGroup.col)),
                              shape=mGroup.shape)
    mGroups.append(mGroup)
    mGroupsNormed.append(mGroupNormed)
    cdGroupNormed = cdGroup / cdGroup.sum()
    cdGroups.append(cdGroup)
    cdGroupsNormed.append(cdGroupNormed)
    # save the matrix
    np.savetxt(foutprefix + ".%d.mat" % (iGroup),
               np.vstack((mGroup.row, mGroup.col, mGroup.data)).T,
               "%6d%6d%6d", header=header+"# bin1 bin2 count\n")
    np.savetxt(foutprefix + ".%d.mat" % (iGroup),
               np.vstack((mGroupNormed.row, mGroupNormed.col, mGroupNormed.data)).T,
               "%6d%6d%6d", header=header+"# bin1 bin2 count\n")
    genomicDist = np.arange(cdGroup.size) * binSize
    np.savetxt(foutprefix + ".%d.contactDecay" % (iGroup),
               np.vstack((genomicDist, cdGroup, cdGroupNormed)).T,
               "%30.15e", header=header+"# GenomicDistance Counts Prob")

# plot contact matrices
nBinsMax = nBinsGroups.max()
plt.clf()
plt.cla()
fig, axes = plt.subplots(math.ceil(nGroups / 5), 5, figsize=(40, 40))
for iGroup in range(nGroups):
    mGroupNormed = mGroupsNormed[iGroup]
    iPlotRow = iGroup // 5
    iPlotCol = iGroup % 5
    axis = axes[iPlotRow, iPlotCol]
    axis.set_aspect('equal')
    im = axis.imshow(mGroupNormed.todense(), interpolation='None', aspect='equal',
                     cmap='plasma',
                     norm=LogNorm(mGroupNormed.data.min(),
                                  max(mGroupNormed.data.min()*10,
                                      mGroupNormed.data.max())))
    cb = fig.colorbar(im, ax=axis)
    cb.set_label("Contact Probability")
    axis.set_xlim(0, nBinsMax)
    axis.set_ylim(0, nBinsMax)
    axis.set_xlabel("Bin (bin size = %5d)" % (binSize))
    axis.set_ylabel("Bin (bin size = %5d)" % (binSize))
    axis.set_title("Gene %5d - %5d" %
                   (iGroup * every, min((iGroup+1) * every, bed.shape[0]) - 1))
for i in range(nGroups, axes.flatten().size):
     fig.delaxes(axes.flatten()[i])
plt.savefig(foutprefix+".mat.png", dpi=300)

# plot contact decay
plt.close(fig)
figureCD = plt.figure()
axCD = figureCD.add_subplot(111)
axCD.set_xlabel("Genomic distance (bp)")
axCD.set_ylabel("Contact probability")
axCD.set_ylim(0, 0.10)
axCD.set_xlim(binSize, 1000)
for iGroup in range(nGroups):
    cdGroupNormed = cdGroupsNormed[iGroup]
    genomicDist = np.arange(cdGroupNormed.size) * binSize
    axCD.plot(genomicDist, cdGroupNormed, label="Gene group %d" % (iGroup),
              c=plt.get_cmap('plasma').colors[int(256 / nGroups * iGroup)],
              linewidth=0.2)
figureCD.savefig(foutprefix+".cd.png", dpi=300)
