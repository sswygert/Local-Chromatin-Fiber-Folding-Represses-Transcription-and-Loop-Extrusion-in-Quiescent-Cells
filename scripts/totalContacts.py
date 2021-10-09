#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2020 dejunlin <dejun.lin@gmail.com>
# Created at 2020-02-27 14:05 on dejunlin@threonine.gs.washington.edu
# Usage: totalContacts.py
# Description: Get total number of contacts from a mcool or cool file
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
from utils import readMcool, cool2pixels

parser = argparse.ArgumentParser()
parser.add_argument("fmcool", type=str,
                    help="Cooler multiple-binsize contact files")
parser.add_argument("--binSize", type=int, default=-1,
                    help="Use this to select the bin size from the input mcool\
                    file. Default to -1, meaning that the inputs are treated as\
                    single-binsize .cool files. For mcool file, the sum of contacts\
                    is invariant to bin size choice.")

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

cool, binSize = readMcool(fmcool, binSize)
p = cool2pixels(cool)
n = p[:].iloc[:, 2].sum()
print(f"{n:15.9e}")
