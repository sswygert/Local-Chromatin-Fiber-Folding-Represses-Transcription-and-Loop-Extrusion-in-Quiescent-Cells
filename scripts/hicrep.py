#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2020 dejunlin <dejun.lin@gmail.com>
# Created at 2020-01-27 15:21 on dejunlin@threonine.gs.washington.edu
# Usage: hicrep.py
# Description: Compute HiCRep reproducibility stratum-corrected correlation score (SCCS). 
# Reference: Genome Res. 2017 Nov;27(11):1939-1949. doi: 10.1101/gr.220640.117
# The algorithm first normalizes the input contact matrices by the total 
# number of contacts and then for each chromosome: 1) mean-filter the input
# matrices with an input window size; 2) exclude common zero entries in
# the input matrices; 3) compute the SCC score. It doesn't have the
# procedure to bootstrap the window-size parameter
#
# Distributed under terms of the MIT license.
import os
import scipy as sp
import numpy as np
import math
import sys
import warnings
import pandas as pd
import cooler
import h5py
import scipy.sparse as sp
import scipy.ndimage as spi
from scipy.stats import rankdata
from statsmodels.distributions.empirical_distribution import ECDF
from utils import readMcool, cool2pixels, getSubCoo, trimDiags

def meanFilterSparse(a: sp.coo_matrix, h: int):
    """Apply a mean filter to an input sparse matrix. This convolves
    the input with a kernel of size 2*h + 1 with constant entries and
    subsequently reshape the output to be of the same shape as input

    Args:
        a: `sp.coo_matrix`, Input matrix to be filtered
        h: `int` half-size of the filter

    Returns:
        `sp.coo_matrix` filterd matrix
    """
    assert h > 0, "meanFilterSparse half-size must be greater than 0"
    assert sp.issparse(a) and a.getformat() == 'coo',\
        "meanFilterSparse input matrix is not scipy.sparse.coo_matrix"
    assert a.shape[0] == a.shape[1],\
        "meanFilterSparse cannot handle non-square matrix"
    fSize = 2 * h + 1
    # filter is a square matrix of constant 1 of shape (fSize, fSize)
    shapeOut = np.array(a.shape) + fSize - 1
    # This is the non-zero blocks of the doubly Toeplitz matrix.
    # For square input matrix a and a kernel of constant, such as the
    # mean filter kernel, mToeplitz @ a @ mToeplitz.T == 
    # doublyToeplitz @ aLinear, where aLinear is linearized a with 
    # (a.row[n-1], a.row[n-2], ..., a.row[0]).T and doublyToeplitz
    # is a doubly Toeplitz matrix formed by H = mToeplitz
    #                  H 0 . . . . 0 
    #                  H H 0 . . . 0 
    #                  . H H 0 . . 0 
    # doublyToeplitz = H . H H 0 . 0 
    #                  0 H . . . 0 0 
    #                  . 0 H . . H 0 
    #                  0 . 0 H . H H 
    # i.e., doublyToeplitz @ aLinear == mToeplitz @ a @ mToeplitz.T 
    # see my previous commit for details
    mToeplitz = sp.diags(np.ones(fSize),
                         np.arange(-fSize+1, 1),
                         shape=(shapeOut[1], a.shape[1]),
                         format='csr')
    ans = sp.coo_matrix((mToeplitz @ a) @ mToeplitz.T)
    # remove the edges since we don't care about them if we are smoothing
    # the matrix itself
    ansNoEdge = ans.tocsr()[h:(h+a.shape[0]), h:(h+a.shape[1])].tocoo()

    # # verify answer here
    # from scipy import signal
    # import scipy.ndimage as spi
    # k = np.ones((fSize, fSize), dtype=float)
    # # compare to full convolution
    # ansLib = signal.convolve2d(a.todense(), k, "full")
    # diff = sp.coo_matrix(ansLib) - ans
    # assert np.isclose(diff.data, np.zeros(diff.nnz, dtype=float)).all()
    # # compare non-edge version
    # ansLibNoEdge = spi.convolve(a.todense(), k, mode='constant', cval=0.0)
    # diffNoEdge = ansNoEdge - sp.csr_matrix(ansLibNoEdge)
    # assert np.isclose(diffNoEdge.data, np.zeros(diffNoEdge.nnz, dtype=float)).all()

    # Assign different number of neighbors to the edge to better
    # match what the original R implementation of HiCRep does
    rowDist2Edge = np.minimum(ansNoEdge.row, ansNoEdge.shape[0] - 1 - ansNoEdge.row)
    nDim1 = h + 1 + np.minimum(rowDist2Edge, h)
    colDist2Edge = np.minimum(ansNoEdge.col, ansNoEdge.shape[1] - 1 - ansNoEdge.col)
    nDim2 = h + 1 + np.minimum(colDist2Edge, h)
    nNeighbors = nDim1 * nDim2
    ansNoEdge.data /= nNeighbors

    # # verify answer here 
    # from scipy import signal
    # import scipy.ndimage as spi
    # kSumNeighbors = np.ones((fSize, fSize), dtype=float)
    # sumNeigbhors = spi.convolve(a.todense(), kSumNeighbors, mode='constant', cval=0.0)
    # nNeighbors = spi.convolve(np.ones_like(sumNeigbhors, dtype=float),
    #                           kSumNeighbors, mode='constant', cval=0.0)
    # ansLibNoEdge = sumNeigbhors / nNeighbors
    # diffNoEdge = ansNoEdge - sp.coo_matrix(ansLibNoEdge)
    # assert np.isclose(diffNoEdge.data, np.zeros(diffNoEdge.nnz, dtype=float)).all()
    return ansNoEdge


def vstran(a: np.ndarray):
    """compute the empirical CDF at each ranked element
    this is the so-called variance stabilizing xformation
    scipy.stats.rankdata doesn't have a method to randomly break
    the original order of tied elements, unlike the R implementation,
    so I use numpy quicksort, which could break the ties somewhat
    randomly

    Args:
        a: `np.ndarray` input array to be transformed

    Returns:
        `np.ndarray` transformed array
    """
    aRank = np.argsort(a, kind='quicksort') + 1
    return ECDF(aRank)(aRank)

def resample(m: sp.coo_matrix, size: int):
    """Resample with replacement the input matrix so that the
    resulting matrix sum to the given size
    Args:
        m: `sp.coo_matrix` Input matrix
        size: Resulting matrix sum to this number

    Returns:
        resampled matrix
    """
    bins = np.arange(m.data.size)
    p = m.data / m.data.sum()
    samples = np.random.choice(bins, size=size, p=p)
    sampledData = np.bincount(samples, minlength=bins.size)
    ans = sp.coo_matrix((sampledData, (m.row, m.col)), shape=m.shape)
    ans.eliminate_zeros()
    return ans

def sccOfDiag(diag1: np.ndarray, diag2: np.ndarray):
    """Get the correlation coefficient and weight of two input
    diagonal arrays

    Args:
        diag1: `np.ndarray` input array 1
        diag2: `np.ndarray` input array 2

    Returns:
       tuple of 2 floats, the Pearson's correlation rho and weight
    """
    # remove common zeros
    idxNZ = np.where((diag1 != 0.0) | (diag2 != 0.0))[0]
    iN = idxNZ.size
    if iN <= 2:
        return (np.nan, np.nan)
    iDiagNZ1 = diag1[idxNZ]
    iDiagNZ2 = diag2[idxNZ]
    # variance stabilizing xformation -- these are for geometric 
    # average of the variances only
    iDiagVS1 = vstran(iDiagNZ1)
    iDiagVS2 = vstran(iDiagNZ2)
    rho = np.corrcoef(iDiagNZ1, iDiagNZ2)[0, 1]
    # the original R implementation impose Bessel's correction
    # in the variance calculation, which I don't think it make sense
    # here
    ws = iN * np.sqrt(np.var(iDiagVS1, ddof=1)*np.var(iDiagVS2, ddof=1))
    return (rho, ws)

def hicrepSCC(cool1: cooler.api.Cooler, cool2: cooler.api.Cooler,
              h: int, dBPMax: int, bDownSample: bool):
    """Compute hicrep score between two input Cooler contact matrices

    Args:
        cool1: `cooler.api.Cooler` Input Cooler contact matrix 1
        cool2: `cooler.api.Cooler` Input Cooler contact matrix 2
        h: `int` Half-size of the mean filter used to smooth the
        input matrics
        dBPMax `int` Only include contacts that are at most this genomic
        distance (bp) away
        bDownSample: `bool` Down sample the input with more contacts
        to the same number of contacts as in the other input

    Returns:
        `float` scc scores for each chromosome
    """
    binSize1 = cool1.info['bin-size']
    binSize2 = cool2.info['bin-size']
    assert binSize1 == binSize2,\
        f"Input cool files {fmcool1} and {fmcool2} have different bin sizes"
    assert cool1.info['nbins'] == cool2.info['nbins'],\
        f"Input cool files {fmcool1} and {fmcool2} have different number of bins"
    assert cool1.info['nchroms'] == cool2.info['nchroms'],\
        f"Input cool files {fmcool1} and {fmcool2} have different number of chromosomes"
    assert (cool1.chroms()[:] == cool2.chroms()[:]).all()[0],\
        f"Input file {fmcool1} and {fmcool2} have different chromosome names"
    binSize = binSize1
    assert dBPMax > binSize, f"Input dBPmax is smaller than binSize"
    # this is the exclusive upper bound
    dMax = dBPMax // binSize + 1
    p1 = cool2pixels(cool1)
    p2 = cool2pixels(cool2)
    bins1 = cool1.bins()
    bins2 = cool2.bins()
    # get the total number of contacts as normalizing constant
    n1 = p1[:]['count'].sum()
    n2 = p2[:]['count'].sum()
    chrNames = cool1.chroms()[:]['name'].to_numpy()
    # filter out mitochondria chromosome
    chrNames = np.array([name for name in chrNames if name != 'M'])
    scc = np.full(chrNames.shape[0], -2.0)
    for iChr in range(chrNames.shape[0]):
        chrName = chrNames[iChr]
        # normalize by total number of contacts
        mS1 = getSubCoo(p1, bins1, chrName)
        assert mS1.size > 0, "Contact matrix 1 of chromosome %s is empty" % (chrName)
        assert mS1.shape[0] == mS1.shape[1],\
            "Contact matrix 1 of chromosome %s is not square" % (chrName)
        mS2 = getSubCoo(p2, bins2, chrName)
        assert mS2.size > 0, "Contact matrix 2 of chromosome %s is empty" % (chrName)
        assert mS2.shape[0] == mS2.shape[1],\
            "Contact matrix 2 of chromosome %s is not square" % (chrName)
        assert mS1.shape == mS2.shape,\
            "Contact matrices of chromosome %s have different input shape" % (chrName)
        nDiags = mS1.shape[0] if dMax < 0 else min(dMax, mS1.shape[0])
        rho = np.full(nDiags, np.nan)
        ws = np.full(nDiags, np.nan)
        # remove major diagonal and all the diagonals >= nDiags
        # to save computation time
        m1 = trimDiags(mS1, nDiags, False)
        m2 = trimDiags(mS2, nDiags, False)
        if bDownSample:
            # do downsampling
            size1 = m1.sum()
            size2 = m2.sum()
            if size1 > size2:
                m1 = resample(m1, size2).astype(float)
            elif size2 > size1:
                m2 = resample(m2, size1).astype(float)
        else:
            # just normalize by total contacts
            m1 = m1.astype(float) / n1
            m2 = m2.astype(float) / n2
        if h > 0:
            # apply smoothing
            m1 = meanFilterSparse(m1, h)
            m2 = meanFilterSparse(m2, h)
        # ignore the main diagonal iD == 0
        for iD in range(1, nDiags):
            iDiag1 = m1.diagonal(iD)
            iDiag2 = m2.diagonal(iD)
            rho[iD], ws[iD] = sccOfDiag(iDiag1, iDiag2)
        wsNan2Zero = np.nan_to_num(ws, copy=True)
        rhoNan2Zero = np.nan_to_num(rho, copy=True)
        scc[iChr] = rhoNan2Zero @ wsNan2Zero / wsNan2Zero.sum()
    return scc


def main(*args):
    import argparse
    import subprocess
    import re

    np.random.seed(10)

    parser = argparse.ArgumentParser()
    parser.add_argument("fmcool1", type=str,
                        help="First cooler multiple-binsize contact files")
    parser.add_argument("fmcool2", type=str,
                        help="Second cooler multiple-binsize contact files")
    parser.add_argument("fout", type=str,
                        help="Output results to this file. Output format would be\
                        one column of scc scores for each chromosome")
    parser.add_argument("--binSize", type=int, default=-1,
                        help="Use this to select the bin size from the input mcool\
                        file. Default to -1, meaning that the inputs are treated as\
                        single-binsize .cool files")
    parser.add_argument("--h", type=int, default=0,
                        help="Smooth the input contact matrices using a 2d mean\
                        filter with window size of 1 + 2 * value")
    parser.add_argument("--dBPMax", type=int, default=-1,
                        help="Only consider contacts at most this number of bp away\
                        from the diagonal. Default to -1, meaning the entire\
                        contact matrix is used")
    parser.add_argument("--bDownSample", action='store_true', default=False,
                        help="Down sample the input with more contact counts to\
                        the the same number of counts as the other input with less\
                        contact counts. If turned off, the input matrices will be\
                        normalized by dividing the counts by their respective total\
                        number of contacts.")

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

    fmcool1 = args.fmcool1
    fmcool2 = args.fmcool2
    fout = args.fout
    binSize = args.binSize
    h = args.h
    dBPMax = args.dBPMax
    bDownSample = args.bDownSample

    cool1, binSize1 = readMcool(fmcool1, binSize)
    cool2, binSize2 = readMcool(fmcool2, binSize)

    scc = hicrepSCC(cool1, cool2, h, dBPMax, bDownSample)

    np.savetxt(fout, scc, "%30.15e", header=header)

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
