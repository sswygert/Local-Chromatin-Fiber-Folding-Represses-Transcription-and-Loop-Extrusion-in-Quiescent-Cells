#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2019 dejunlin <dejun.lin@gmail.com>
# Created at 2019-10-14 12:07 on dejunlin@threonine.gs.washington.edu
# Usage: selectBins.py
# Description: Get a list of bin ids for the input list of genomic loci
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

parser = argparse.ArgumentParser()
parser.add_argument("finput", type=str,
                    help="Input a list loci in the format: label chr start end")
parser.add_argument("fmcool", type=str,
                    help="Input cooler multiple-binsize contact files\
                    whose bins are used for the selection")
parser.add_argument("binSize", type=int,
                    help="Use data corresponding to this bin size in fmcool\
                    one of 10, 50, 100, 200, 500, 1000, 3200 and 5000 bp")
parser.add_argument("fout", type=str, default="",
                    help="Output bin index to this file")

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

finput = args.finput
fmcool = args.fmcool
binSize = args.binSize
fout = args.fout

loci = np.loadtxt(finput, dtype=str)

"""For example:
loci = np.array([
    ["LTLMI",   "I",     "1",   "801",  ],
    ["LTLMII",  "II",    "1",   "6608", ],
    ["LTLMIII", "III",   "1",   "1098", ],
    ["LTLMIV",  "IV",    "1",   "904",  ],
    ["LTLMV",   "V",     "1",   "6473", ],
    ["LTLMVI",  "VI",    "1",   "5530", ],
    ["LTLMVII", "VII",   "1",   "781",  ],
    ["LTLMVIII","VIII",  "1",   "5505", ],
    ["LTLMIX",  "IX",    "1",   "7784", ],
    ["LTLMX",   "X",     "1",   "7767", ],
    ["LTLMXI",  "XI",    "1",   "807",  ],
    ["LTLMXII", "XII",   "1",   "12085",],
    ["LTLMXIII","XIII",  "1",   "6344", ],
    ["LTLMXIV", "XIV",   "1",   "7428", ],
    ["LTLMXV",  "XV",    "1",   "847",  ],
    ["LTLMXVI", "XVI",   "1",   "7223", ],
    ["CENMI",   "I",    "151465",   "151582"],
    ["CENMII",  "II",   "238207",   "238323"],
    ["CENMIII", "III",  "114385",   "114501"],
    ["CENMIV",  "IV",   "449711",   "449821"],
    ["CENMV",   "V",    "151987",   "152104"],
    ["CENMVI",  "VI",   "148510",   "148627"],
    ["CENMVII", "VII",  "496920",   "497038"],
    ["CENMVIII","VIII", "105586",   "105703"],
    ["CENMIX",  "IX",   "355629",   "355745"],
    ["CENMX",   "X",    "436307",   "436425"],
    ["CENMXI",  "XI",   "440129",   "440246"],
    ["CENMXII", "XII",  "150828",   "150947"],
    ["CENMXIII","XIII", "268031",   "268149"],
    ["CENMXIV", "XIV",  "628758",   "628875"],
    ["CENMXV",  "XV",   "326584",   "326702"],
    ["CENMXVI", "XVI",  "555957",   "556073"],
    ["RTLMI",   "I",     "229411",   "230218" ],
    ["RTLMII",  "II",    "812379",   "813184" ],
    ["RTLMIII", "III",   "315783",   "316620" ],
    ["RTLMIV",  "IV",   "1524625",   "1531933"],
    ["RTLMV",   "V",     "569599",   "576874" ],
    ["RTLMVI",  "VI",    "269731",   "270161" ],
    ["RTLMVII", "VII",  "1083635",   "1090940"],
    ["RTLMVIII","VIII",  "556105",   "562643" ],
    ["RTLMIX",  "IX",    "439068",   "439888" ],
    ["RTLMX",   "X",     "744902",   "745751" ],
    ["RTLMXI",  "XI",    "665904",   "666816" ],
    ["RTLMXII", "XII",  "1064281",   "1078177"],
    ["RTLMXIII","XIII",  "923541",   "924431" ],
    ["RTLMXIV", "XIV",   "783278",   "784333" ],
    ["RTLMXV",  "XV",   "1083922",   "1091291"],
    ["RTLMXVI", "XVI",   "942396",   "948066" ]
    ])
"""

# Read the mcool file
mcool = h5py.File(fmcool, 'r')
cool = cooler.Cooler(mcool['resolutions'][str(binSize)])

fhout = open(fout, 'w')
fhout.write(header)
for locus in loci:
    region = locus[1] + ":" + locus[2] + "-" + locus[3]
    mBins = cool.bins().fetch(region)
    bins = mBins.index.values
    fhout.write("%10s" % (locus[0]))
    for iBin in bins:
        fhout.write("%10d" % iBin)
    fhout.write("\n")
fhout.close()
