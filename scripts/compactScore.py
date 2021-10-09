#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2019 dejun <dejun.lin@gmail.com>
# Created at 2019-08-19 20:42 on dejun@dejun-GS60-2QE
# Usage: compactScore.py
# Description: Compute the compact score for all the budding yeast genes as 
# shown in Fig. S6 of https://doi.org/10.1016/j.cell.2015.05.048
#
# Distributed under terms of the MIT license.
import os
import scipy as sp
from scipy.stats import norm
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
from intermine.webservice import Service
import gffutils
import swifter
from timeit import default_timer as timer
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator, LogFormatterSciNotation
from sklearn.neighbors import KDTree


parser = argparse.ArgumentParser()
parser.add_argument("fmcool", type=str,
                    help="Input cooler multiple-binsize contact files")
parser.add_argument("binSize", type=int,
                    help="Use data corresponding to this bin size in fmcool\
                    one of 200, 500, 1000, 3200 and 5000 bp")
parser.add_argument("fout", type=str,
                    help="Output results to this file: geneSymbol geneName chr start end\
                    totalContacts contactsBelowDiagOne length totalContactsExpected\
                    totalContactsNormed")
parser.add_argument("--fpng", type=str, default="",
                    help="Optionally generate png plotting Fig. S6 A-D")
parser.add_argument("--fgtf", type=str, default="",
                    help="Compute the score only for the ORF in this gtf file instead of\
                    all genes from YestMine")

args = parser.parse_args()

header = "#"+" ".join(sys.argv)+"\n"

gitrev = subprocess.Popen('cd '+os.path.dirname(os.path.realpath(__file__))+' && git rev-parse HEAD --abbrev-ref HEAD', shell=True, stdout=subprocess.PIPE,\
                       stderr=subprocess.PIPE, encoding='utf-8')
if not gitrev.stderr.read():
    p = re.compile("\n")
    header += "# @rev " + p.sub("\n# @branch ", gitrev.stdout.read().strip()) + "\n"

fmcool = args.fmcool
binSize = args.binSize
fout = args.fout
fpng = args.fpng
fgtf = args.fgtf

# Read the mcool file
mcool = h5py.File(fmcool, 'r')
cool = cooler.Cooler(mcool['resolutions'][str(binSize)])

geneList = None
colnames = ['symbol', 'name', 'chr', 'start', 'end']
if fgtf == "":
    # Generate query for all the genes
    service = Service("https://yeastmine.yeastgenome.org:443/yeastmine/service")
    query = service.new_query("Gene")
    query.add_view(
            "primaryIdentifier", "qualifier", "symbol", "name",
            "chromosomeLocation.start", "chromosomeLocation.end",
            "chromosome.primaryIdentifier"
    )
    query.add_sort_order("Gene.chromosome.primaryIdentifier", "ASC")
    query.add_constraint("status", "=", "Active", code = "E")
    query.add_constraint("status", "IS NULL", code = "F")
    query.add_constraint("featureType", "=", "ORF", code = "A")
    query.set_logic("A and (E or F) and E and E")

    # only keep verified genes and geneome chromosome (not mitochrondria)
    # remove 'chr' from chromosome name to be consistent with cooler
    geneList = pd.concat([pd.DataFrame(
        [[row["symbol"], row["name"],
          re.sub(r"chr", "", row["chromosome.primaryIdentifier"]),
          row["chromosomeLocation.start"], row["chromosomeLocation.end"]]])
        for row in query.rows()
        if row["qualifier"] == "Verified" and
        row["chromosome.primaryIdentifier"] != "chrmt"] )

else:
    fdb = fgtf+'.db'
    db = gffutils.create_db(fgtf, dbfn=fdb, force=True, keep_order=True,
                            merge_strategy='merge', sort_attribute_values=True,
                            disable_infer_transcripts=True)
    geneList = pd.concat([pd.DataFrame(
        [["N/A", gene.id, re.sub(r"chr", "", gene.chrom), gene.start, gene.end]])
        for gene in db.all_features()
        ])
if geneList is None:
    raise ValueError("Can't not create list of genes and their genomic coordinates")
geneList.columns = colnames
geneList.sort_values(by=['chr','start','end'],inplace=True)
#Remove spaces from name 
geneList['name'] = geneList['name'].str.title().replace(r"\s+", r"", regex=True)

#sum the contact counts in the body of each gene
#Total counts
geneList['totalContacts'] = geneList.apply(
    lambda gene: cool.matrix(as_pixels=True, balance=False, sparse=True).fetch(
        gene['chr'] + ':' + str(gene['start']) + '-' + str(gene['end']))['count'].sum(),
        axis=1)
#Below a certain diagonal
def sumBelowDiagonal(gene, maxDiag):
    m = cool.matrix(as_pixels=True, balance=False, sparse=True).fetch(
        gene['chr'] + ':' + str(gene['start']) + '-' + str(gene['end']))
    return m.loc[m['bin2_id'] - m['bin1_id'] <= maxDiag]['count'].sum()
geneList['contactsBelowDiagOne'] = geneList.apply(sumBelowDiagonal, maxDiag=1, axis=1)

#compute the gene length
geneList['length'] = geneList['end'] - geneList['start'] + 1

#compute the k-nn for normalization
f = geneList.loc[:, ['contactsBelowDiagOne','length']].to_numpy()
#min-max normalize 
fMins = f.min(axis=0)
fMaxs = f.max(axis=0)
fNorm = (f - f.min(axis=0)) / (fMaxs - fMins)
v = geneList.loc[:, ['totalContacts']].to_numpy()
kdTree = KDTree(fNorm)
#TODO: exclude the query point itself because it is part of the input
nearest_dist, nearest_ind = kdTree.query(fNorm, k=20+1)
#compute an gaussian-weighted average 
normVar = nearest_dist.var(axis=1)
weights = np.zeros(nearest_dist.shape)
for i in range(weights.shape[0]):
    weights[i] = norm.pdf(nearest_dist[i], loc=0, scale=1.0)

geneList['totalContactsExpected'] = np.average(v[nearest_ind][:,:,0], axis=1, weights=weights)
geneList['totalContactsNormed'] = geneList['totalContacts'] / geneList['totalContactsExpected']

# save to file -- this can be read back in bypd.read_table('filename.txt', sep=r'\s+') 
with open(fout,'w') as outfile:
    outfile.write(header)
    geneList.to_string(outfile, float_format="%30.15e")

locator = LogLocator(base=2)
formatter = LogFormatterSciNotation(base=2)
def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax, ticks=locator, format=formatter)

if fpng != "":
    mpl.rcParams.update({'font.size': 25})
    fig, ax = plt.subplots(2, 2, figsize=(20,20))
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.3)
    fig.suptitle('Compactness score for all genes in S. cerevisiae')
    ###
    plt.cla()
    # unormed
    ax[0, 0].set_xlabel('# of contacts below 2nd diagonal')
    ax[0, 0].set_ylabel('Gene Length')
    ax[0, 0].set_xscale('log', basex=2)
    ax[0, 0].set_yscale('log', basey=2)
    ax[0, 0].set_xlim(1, geneList['contactsBelowDiagOne'].max())
    ax[0, 0].set_ylim(geneList['length'].min(), geneList['length'].max())
    sc0 = ax[0, 0].scatter(geneList['contactsBelowDiagOne'], geneList['length'],
                    c=geneList['totalContacts'], s=1, cmap='inferno',
                    norm=LogNorm())
    cb0 = colorbar(sc0)
    cb0.set_label('# contacts within gene')
    # expected
    ax[0, 1].set_xlabel('# of contacts below 2nd diagonal')
    ax[0, 1].set_ylabel('Gene Length')
    ax[0, 1].set_xscale('log', basex=2)
    ax[0, 1].set_yscale('log', basey=2)
    ax[0, 1].set_xlim(1, geneList['contactsBelowDiagOne'].max())
    ax[0, 1].set_ylim(geneList['length'].min(), geneList['length'].max())
    sc1 = ax[0, 1].scatter(geneList['contactsBelowDiagOne'], geneList['length'],
                    c=geneList['totalContactsExpected'], s=1, cmap='inferno',
                    norm=LogNorm())
    cb1 = colorbar(sc1)
    cb1.set_label('Expected # contacts within gene')

    # normed
    ax[1, 0].set_xlabel('Gene Length')
    ax[1, 0].set_ylabel('Normalized # contacts within gene')
    sc2 = ax[1, 0].scatter(geneList['length'], geneList['totalContactsNormed'], s=1)

    ax[1, 1].set_xlabel('# of contacts below 2nd diagonal')
    ax[1, 1].set_ylabel('Normalized # contacts within gene')
    sc3 = ax[1, 1].scatter(geneList['contactsBelowDiagOne'], geneList['totalContactsNormed'], s=1)
    ###
    fig.savefig(fpng, dpi=300)
