#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2019 dejunlin <dejun.lin@gmail.com>
# Created at 2019-09-12 11:03 on dejunlin@threonine.gs.washington.edu
# Usage: mcool2coo.py
# Description: Convert a multi-binsize hi-c mcool file to a coo matrix
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
import cooler
import h5py
from scipy.sparse import coo_matrix

parser = argparse.ArgumentParser()
parser.add_argument("fmcool", type=str,
                    help="Input mcool file")
parser.add_argument("binSize", type=int,
                    help="The binSize in the fmcool to output")
parser.add_argument("fcoo", type=str,
                    help="Output coo file")
parser.add_argument("--bBalance", action='store_true', default=False,
                    help="If apply matrix balancing (or ICE normalization) to output")

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
fcoo = args.fcoo
bBalance = args.bBalance

mcool = h5py.File(fmcool, 'r')
cool = cooler.Cooler(mcool['resolutions'][str(binSize)])

coo = cool.matrix(balance=False, sparse=True)[:,:]
if bBalance:
    # the parameters were taken from Seungsoo's pipeline as reported on 
    # https://doi.org/10.1016/j.molcel.2018.11.020
    bias, stats = cooler.balance_cooler(cool, mad_max=9, ignore_diags=1)
    # Cooler fill in NaN for bins that are filtered out according to
    # https://github.com/mirnylab/cooler/blob/843dadca5ef58e3b794dbaf23430082c9a634532/cooler/balance.py#L333
    # To prevent propagating NaN downstream, I replace them with zero here
    np.nan_to_num(bias, copy=False)
    coo.data = bias[coo.row] * bias[coo.col] * coo.data

np.savetxt(fcoo, np.vstack((coo.row, coo.col,coo.data)).T, "%10d%10d%30.15e", header=header)
