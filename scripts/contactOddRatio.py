#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2021 dejun <dejun.lin@gmail.com>
# Created at 2021-09-23 18:35 on dejun@Z590E
# Usage: contactOddRatio.py
# Description: Compute the ratio of contacts
# probability below and above a genomic distance
# threshold
#
# Distributed under terms of the MIT license.
import os
import numpy as np
import sys
import argparse
import subprocess
import re

parser = argparse.ArgumentParser()
parser.add_argument("finput", type=str,
                    help="Input contact probability file. The first column is \
                    the genomic distance and the 3rd column is the contact \
                    probability")
parser.add_argument("foutput", type=str,
                    help="Output the probabilities and odd ratio to this file")
parser.add_argument("min1", type=int,
                    help="Genomic distance lower bound of the denominator\
                    contact probability")
parser.add_argument("max1", type=int,
                    help="Genomic distance upper bound of the denominator\
                    contact probability")
parser.add_argument("min2", type=int,
                    help="Genomic distance lower bound of the numerator\
                    contact probability")
parser.add_argument("max2", type=int,
                    help="Genomic distance upper bound of the numerator\
                    contact probability")

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
foutput = args.foutput
min1 = args.min1
max1 = args.max1
min2 = args.min2
max2 = args.max2

probs = np.loadtxt(finput, usecols=(0, 2))
bins1 = np.where(np.logical_and(min1 <= probs[:, 0], probs[:, 0] <= max1))[0]
bins2 = np.where(np.logical_and(min2 <= probs[:, 0], probs[:, 0] <= max2))[0]

p1 = probs[bins1, 1].sum()
p2 = probs[bins2, 1].sum()
r = p2 / p1
out = np.array([[p1, p2, r]])

header = header + "# prob1 prob2 prob2/prob1"
np.savetxt(foutput, out, "%30.15e", header=header)
