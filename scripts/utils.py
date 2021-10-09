#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2019 dejunlin <dejun.lin@gmail.com>
# Created at 2019-11-14 10:09 on dejunlin@threonine.gs.washington.edu
# Usage: utils.py
# Description: utility functions and classes
#
# Distributed under terms of the MIT license.
import numpy as np
import pandas as pd
import cooler
import h5py
import math
import scipy.sparse as sp
from scipy.sparse import coo_matrix

def readMcool(fmcool: str, binSize: int):
    """Read from a mcool or cool file and return the Cooler object

    Args:
        fmcool: Input file name
        binSize: Bin size to select from the mcool file. If this value
        is <= 0, the input will be treated as a cool file instead

    Returns:
        cooler.api.Cooler object
    """
    mcool = h5py.File(fmcool, 'r')
    if binSize > 0:
        return cooler.Cooler(mcool['resolutions'][str(binSize)]), binSize
    else:
        cool = cooler.Cooler(mcool)
        return cool, cool.binsize

def pixels2Coo(df: pd.DataFrame, bins: pd.DataFrame):
    """Convert Cooler's contact matrix in "pixels" DataFrame to
    scipy coo_matrix. The "pixels" format is a 3-column DataFrame:
    'bin1_id', 'bin2_id', 'counts' for each unique contact

    Args:
        df: Input DataFrame
        bins: Cooler bins for the contacts

    Returns:
        coo_matrix of the input
    """
    binOffset = bins.index[0]
    nBins = bins.shape[0]
    df['bin1_id'] -= binOffset
    df['bin2_id'] -= binOffset
    return coo_matrix((df['count'].to_numpy(), (df['bin1_id'].to_numpy(),
                                                df['bin2_id'].to_numpy())),
                      shape=(nBins, nBins))

def cool2pixels(cool: cooler.api.Cooler):
    """Return the contact matrix in "pixels" format

    Args:
        cool: Input cooler object

    Returns:
        cooler.core.RangeSelector2D object
    """
    return cool.matrix(as_pixels=True, balance=False, sparse=True)

def trimDiags(a: sp.coo_matrix, iDiagMax: int, bKeepMain: bool):
    """Remove diagonal elements whose diagonal index is >= iDiagMax
    or is == 0

    Args:
        a: Input scipy coo_matrix
        iDiagMax: Diagonal offset cutoff
        bKeepMain: If true, keep the elements in the main diagonal;
        otherwise remove them

    Returns:
        coo_matrix with the specified diagonals removed
    """
    gDist = np.abs(a.row - a.col)
    idx = np.where((gDist < iDiagMax) & (bKeepMain | (gDist != 0)))
    return sp.coo_matrix((a.data[idx], (a.row[idx], a.col[idx])),
                         shape=a.shape, dtype=a.dtype)


def getSubCoo(pixels: cooler.core.RangeSelector2D, bins: cooler.core.RangeSelector1D,
              regionStr: str):
    """Fetch a region from Cooler a contact matrix and return it as a
    coo_matrix

    Args:
        pixels: Input Cooler range selector object of the contact matrix
        bins: Input Cooler range selector object of the bin definition
        regionStr: String for selecting genomic region

    Returns:
        coo_matrix contact matrix corresponding to the input region
    """
    mSub = pixels.fetch(regionStr)
    # Assume Cooler always use upper triangle
    assert (mSub['bin1_id'] <= mSub['bin2_id']).all(),\
        f"Contact matrix of region {regionStr} has lower-triangle entries"
    binsSub = bins.fetch(regionStr)
    return pixels2Coo(mSub, binsSub)

def interpolateArray(a, n):
    """Interpolate the elements of an array so that its length
    is n

    Args:
        a: `np.ndarray` input array
        n: `int` length of output array

    Returns:
    """
    xInput = np.linspace(0, 1.0, a.shape[0])
    xOutput = np.linspace(0, 1.0, n)
    return np.interp(xOutput, xInput, a)

"""
Compute the contact decay for a genomic region specified by 'region',
which is something like: 'I:1-100', the first 100 bp of chromosome I
"""
def contactDecay(cool, region, normLength: int = 0, bNormNpairs: bool = False,
                 nucPeaks: pd.DataFrame = None):
    """Compute total number of contacts at each genomic distance across an
    input region

    Args:
        cool: input Cooler handle of contact matrices
        region: `string` input genomic region of interest of the form 'chrom:start-end'
        normLength: `int` when none zero, normalize the output array to be
        of length of this number
        bNormNpairs: `bool` whether to normalize the contact counts at each genomic
        distance by the total number of genomic loci pairs at that distance
        nucPeaks: `pd.DataFrame` if not None, compute the inter-nucleosomal contact
        decay with nucleosome locations in this DataFrame: chr loc
        where loc is 0-based location in bp of the nucleosome

    Returns:
        `np.ndarray` total number of contacts at each possible genomic distance
        bin or each nucleosome bin if nucPeaks is specified. In the genomic distance
        bin case, the bin size is the same as in the input cool handle and the
        total number of bins equals to the number of bins in the input region.
        This is aligned such that the 1st element is always the intra-bin contact,
        the 2nd element is always the +1 bin contacts and etc. but different input
        region can give different number of elements in the output unless normLength
        is specified where different region always give normLength-size output
    """
    m = cool.matrix(as_pixels=True, balance=False, sparse=True).fetch(region)
    binsCooler = cool.bins().fetch(region)
    nBinsTotal = binsCooler.shape[0]
    if not m.size:
        return np.empty(0, dtype=float)
    if nucPeaks is not None:
        # convert genomic position-based bins to nucleosome bins
        locMax = binsCooler['end'].tail(1).iloc[0]
        nucPeaksSub = nucPeaks[nucPeaks['chr'] == 'chr' + region].copy().reset_index(drop=True)
        # allow +/- 100 bp around the peak to be called a nucleosome
        nucPeaksSub['locMin'] = (nucPeaksSub['loc'] - 100).clip(lower=0)
        nucPeaksSub['locMax'] = (nucPeaksSub['loc'] + 100).clip(upper=locMax-1)
        # for overlapping nucleasome peak regions, take the average between
        # the two immediate neighboring peaks
        iOverlapLeft = (nucPeaksSub['locMax'].shift(1) >= nucPeaksSub['locMin'])
        iOverlapRight = iOverlapLeft.shift(-1, fill_value=False)
        locOverlapMean = (nucPeaksSub.loc[iOverlapLeft, 'locMin'].reset_index(drop=True) +
         nucPeaksSub.loc[iOverlapRight, 'locMax'].reset_index(drop=True)) // 2
        locOverlapMean.index = np.where(iOverlapLeft)[0]
        nucPeaksSub.loc[iOverlapLeft, 'locMin'] = locOverlapMean
        locOverlapMean.index = np.where(iOverlapRight)[0]
        nucPeaksSub.loc[iOverlapRight, 'locMax'] = locOverlapMean
        # cooler's bin definition is left side inclusive and right side exclusive
        bins = pd.IntervalIndex.from_arrays(binsCooler['start'], binsCooler['end'], closed='left')
        # assign each nucleosome region to a range of Cooler bins
        nucPeaksSub['binMin'] = pd.cut(nucPeaksSub['locMin'], bins, labels=False).\
            apply(lambda df : bins.get_loc(df))
        nucPeaksSub['binMax'] = pd.cut(nucPeaksSub['locMax'], bins, labels=False).\
            apply(lambda df : bins.get_loc(df))
        # for each contact pair in Cooler, assign a nucleosome pair
        binsNuc = pd.IntervalIndex.from_arrays(nucPeaksSub['binMin'], nucPeaksSub['binMax'],
                                               closed='left')
        m['bin1_id'] = pd.cut(m['bin1_id']-binsCooler.index[0], binsNuc, labels=False).\
            apply(lambda df : binsNuc.get_loc(df))
        m['bin2_id'] = pd.cut(m['bin2_id']-binsCooler.index[0], binsNuc, labels=False).\
            apply(lambda df : binsNuc.get_loc(df))
        # drop contact pairs that aren't assigned to a nucleosome
        m = m.loc[~pd.isna(m['bin1_id']) &
                  ~pd.isna(m['bin2_id'])].reset_index(drop=True)
        m['bin1_id'] = m['bin1_id'].astype(int)
        m['bin2_id'] = m['bin2_id'].astype(int)
        # this is the total number of bins possible
        nBinsTotal = nucPeaksSub.shape[0]

    """
    compute the bin distance for each contact
    TODO: this might not be exactly the contacts within genes because
    binning countacts leaves certain read counts outside the gene body
    """
    m['binDist'] = m.apply(lambda c: c.bin2_id - c.bin1_id, axis=1)
    # group by bin distance and sum countacts in each group
    # which results in the histogram as a function of number of bins: 0, 1, ...
    hist = m.groupby('binDist').sum()
    cd = hist['count'].to_numpy()
    binDists = hist.index.to_numpy()
    cdWhole = np.zeros(nBinsTotal, dtype=float)
    cdWhole[binDists] = cd
    if bNormNpairs:
        nPairs = nBinsTotal - binDists
        cdWhole[binDists] = cdWhole[binDists] / nPairs
    if normLength:
        return interpolateArray(cdWhole, normLength)
    return cdWhole

def chipTrack(bigWig, gene, binSize=None):
    """Get the ChIP signal track for a gene

    Args:
        bigWig: Input pyBigWig object as handle of the ChIP BigWig data
        gene: Target gene
        binSize: bin size of the bigWig in bp. If omitted, use the size
        of the 1st bin in bigWig

    Returns:
       `np.ndarray` the 1d ChIP track of the target gene
    """
    if binSize is None:
        # if no input binsize, use the size of the first bin
        bw1stBin = bigWig.intervals(list(bigWig.chroms().keys())[0])[0]
        bs = bw1stBin[1] - bw1stBin[0]
    else:
        bs = binSize
    chip = bigWig.intervals(gene['chrom'], gene['start'], gene['end'])
    iChipTrack = []
    lastBin = 0
    nextBin = 99999999
    for iChip in range(len(chip)):
        iChipStart = chip[iChip][0]
        iChipEnd = chip[iChip][1]
        if (iChipEnd - iChipStart) % bs:
            raise ValueError("BigWig interval %s %d %d is not of multiple of size %d bp" %
                             (gene['chrom'], iChipStart, iChipEnd, bs))
        """Sometimes the BigWig file contains bad bins which are much larger
        than the input binsize but I interpolate them anyway. But I need to
        be careful because this bin can overlap with previous bins or even
        step outside the boundary of the gene"""
        if iChip > 0:
            lastBin = chip[iChip - 1][1]
        else:
            lastBin = iChipStart
        if iChip < len(chip) - 1:
            nextBin = chip[iChip + 1][0]
        else:
            nextBin = iChipEnd
        nBinsChip = int((max(lastBin, iChipEnd) - min(nextBin, iChipStart)) / bs)
        for iBin in range(nBinsChip):
            """Assume all the bins within this pseudo-bin share the same value"""
            iChipTrack.append(chip[iChip][2])
    return np.array(iChipTrack)
