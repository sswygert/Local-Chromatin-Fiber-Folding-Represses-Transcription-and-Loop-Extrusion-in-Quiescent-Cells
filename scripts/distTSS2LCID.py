#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2020 dejun <dejun.lin@gmail.com>
# Created at 2020-10-22 13:06 on dejun@dejun-GS60-2QE
# Usage: distTSS2CID.py
# Description: Get distance of a list of TSS to their
# closest L-CID boundaries
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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument("fTSS", type=str,
                    help="Input bed file for the list of TSS")
parser.add_argument("fCID", type=str,
                    help="Input bed file for the list of L-CID boundaries")
parser.add_argument("fout", type=str,
                    help="Output the distance for each entry in fTSS to the\
                    closest fCID entry")

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

fTSS = args.fTSS
fCID = args.fCID
fout = args.fout

tss = BedTool(fTSS).sort()
cid = BedTool(fCID).sort()

colNames = [ 'col' + str(i) for i in range(23) ]
colNames.append('dist')
out = tss.closest(cid, d=True, t='first')
out.saveas(fout)

