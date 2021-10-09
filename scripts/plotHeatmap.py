#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2020 dejun <dejun.lin@gmail.com>
# Created at 2020-10-14 10:37 on dejun@dejun-GS60-2QE
# Usage: plotHeatmap.py
# Description: Plot heatmap using seaborn.heatmap
#
# Distributed under terms of the MIT license.
import os
import scipy as sp
import scipy.sparse as sps
import numpy as np
import math
import sys
import argparse
import warnings
import subprocess
import re
import matplotlib
import matplotlib.colors as colors
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument("finput", type=str,
                    help="Input matrix file in dense format to be plotted as heatmap")
parser.add_argument("fout", type=str,
                    help="Output heatmap filename")
parser.add_argument("--isCOO", action='store_true', default=False,
                    help="If set, treat the input as COO matrix")
parser.add_argument("--doSymCOO", action='store_true', default=False,
                    help="If symmetrize the input matrix, only do so in case\
                    the input matrix is COO format. Setting this implies isCOO")
parser.add_argument("--tickEvery", type=int, default=1,
                    help="Only show every this number of ticks")
parser.add_argument("--xLabel", type=str, default="Bins",
                    help="x-axis label")
parser.add_argument("--yLabel", type=str, default="Bins",
                    help="y-axis label")
parser.add_argument("--bTicks0OriginCenter", action='store_true', default=False,
                    help="Use the central bin along axis 0 of the input matrix\
                    as origin, i.e., tics would be -nBins//2, ..., 0, ...\
                    nBins//2")
parser.add_argument("--bTicks1OriginCenter", action='store_true', default=False,
                    help="Use the central bin along axis 1 of the input matrix\
                    as origin, i.e., tics would be -nBins//2, ..., 0, ...\
                    nBins//2")
parser.add_argument("--tickLabelMultiplier", type=float, default=1.0,
                    help="Multiply this value to the bin index of the input matrix\
                    to get the x- and y-tick labels. This is only applied to the\
                    axis if bTicks0OriginCenter or bTicks1OriginCenter is set")
parser.add_argument("--cBounds", nargs=2, type=float, default=None,
                    help="Min and max for the output color bar")
parser.add_argument("--cLog", action='store_true', default=False,
                    help="If set, use log scale in the output color bar")
parser.add_argument("--title", type=str, default="",
                    help="Title of the plot")
parser.add_argument("--cLabel", type=str, default="",
                    help="Colorbar label of the plot")
parser.add_argument("--cMap", type=str, default="plasma",
                    help="Matplotlib's cmap for the color scale")
parser.add_argument("--fontSize", type=int, default=20,
                    help="Font size of the plot")

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

finput = args.finput
fout = args.fout
isCOO = args.isCOO
doSymCOO = args.doSymCOO
tickLabelMultiplier = args.tickLabelMultiplier
tickEvery = args.tickEvery
xLabel = args.xLabel
yLabel = args.yLabel
bTicks0OriginCenter = args.bTicks0OriginCenter
bTicks1OriginCenter = args.bTicks1OriginCenter
cBounds = args.cBounds
cLog = args.cLog
title = args.title
cLabel = args.cLabel
cMap = args.cMap
fontSize = args.fontSize

m = np.loadtxt(finput, dtype=float)

if doSymCOO:
    isCOO = True

if isCOO:
    nRows = m[:, :2].astype(int).max() + 1
    m = sps.coo_matrix((m[:, 2].astype(float),
                        (m[:, 0].astype(int), m[:, 1].astype(int))),
                       shape=(nRows, nRows)).toarray()

if doSymCOO:
    m = m.T + m

plt.rcParams.update({'font.size': fontSize})
fig = plt.figure(figsize=(13,10),dpi=200)
ax = plt.gca()
if bTicks0OriginCenter:
    yticks = ((np.arange(m.shape[0]) - m.shape[0] // 2) * tickLabelMultiplier)[tickEvery::tickEvery]
else:
    yticks = (np.arange(m.shape[0]) * tickLabelMultiplier)[tickEvery::tickEvery]

if bTicks1OriginCenter:
    xticks = ((np.arange(m.shape[1]) - m.shape[1] // 2) * tickLabelMultiplier)[tickEvery::tickEvery]
else:
    xticks = (np.arange(m.shape[1]) * tickLabelMultiplier)[tickEvery::tickEvery]

if cBounds is not None:
    vMin = cBounds[0]
    vMax = cBounds[1]
else:
    vMin = m.min()
    vMax = m.max()

if cLog:
    norm = colors.SymLogNorm(linthresh=1e-3, vmin=vMin, vmax=vMax, base=10)
else:
    norm = None

yTickPositions = np.arange(m.shape[0])[tickEvery::tickEvery]
xTickPositions = np.arange(m.shape[1])[tickEvery::tickEvery]
sns.heatmap(m,
            vmin=vMin, vmax=vMax, norm=norm,
            cmap=cMap, cbar_kws={'label' : cLabel},
            ax=ax, linewidths=0.0, rasterized=True
            )

ax.set_xticks(ticks=xTickPositions)
ax.set_yticks(ticks=yTickPositions)
ax.set_xticklabels(labels=xticks, rotation=90)
ax.set_yticklabels(labels=yticks, rotation=0)
ax.set_xlabel(xLabel)
ax.set_ylabel(yLabel)
ax.set_title(title)

fig.savefig(fout, dpi=300, metadata={"Title" : header})
