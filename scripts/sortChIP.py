#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2019 dejunlin <dejun.lin@gmail.com>
# Created at 2019-11-13 10:09 on dejunlin@threonine.gs.washington.edu
# Usage: sortChIP.py
# Description: Sort the genes by the total ChIP occupancy per bp of
# some protein
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
import pyBigWig
import matplotlib
matplotlib.use('Agg')
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10})
from utils import chipTrack

parser = argparse.ArgumentParser()
parser.add_argument("fbw", type=str,
                    help="Input BigWig for ChIP track")
parser.add_argument("fbed", type=str,
                    help="Sort genes in this annotation")
parser.add_argument("foutprefix", type=str,
                    help="Prefix for outputting files")

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

fbw = args.fbw
fbed = args.fbed
foutprefix = args.foutprefix

# Read the bed file
bed = BedTool(fbed).to_dataframe()
# remove '#' from the 1st column name
bed.columns = bed.columns.str.replace('.*chrom$', 'chrom', regex=True)
# Read the BigWig file
bw = pyBigWig.open(fbw)

# go through each gene and sum chip
sumChIPs = np.zeros(bed.shape[0], dtype=float)
for index, gene in bed.iterrows():
    geneLength = gene['end'] - gene['start']
    sumChIPs[index] = chipTrack(bw, gene).sum() / geneLength
bed['ChIP'] = sumChIPs
# sort the chip
colnames = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart',
            'thickEnd', 'itemRGB', 'blockCount', 'blockSizes', 'blockStart', 'ChIP']
sortedBed = bed.sort_values(by=['ChIP'], ascending=False,
                            kind='mergesort'
                            ).loc[:, colnames]
colnames[0] = "#chrom"
BedTool.from_dataframe(sortedBed).saveas(foutprefix+'.bed',
                                      trackline=str.join("\t", colnames))
