#!/usr/bin/env nextflow
nextflow.preview.dsl=2

params.inputCools = "./micro-c/{Q_from_paper,Q_R17R19A_merged,Q_5toA_merged,Q_TSA_merged}/*.NoIN.pair.*.mcool"
params.inputPeaks = "./LCID/Q*MACS1.5.bed"
params.outdir = "./LCID_micro-c_pileup"
params.binSize = 200
params.nBins = 100
params.gDistNBins = 40000 / params.binSize
params.bCmpPeaks2Ways = false
// if we want to match the cool and peak files by sample names
params.bMatchCoolsAndPeaks = false
params.cvmin = 0
params.cvmax = 1
params.bClog = false

/* include './utils' */
include { getPrefix; getMicroCSampleName } from './utils'

process pileupContacts {
  executor = 'slurm'
  queue = 'campus-new'
  memory = '22 GB'
  cpus = 8

  cache 'lenient'

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  tuple file(fmcool), file(fbed), val(excludeZeroSub), val(bMedian), val(bNormDiag), val(bIntra), val(bInter)

  output:
  tuple file(outCoo), file(outSVG), emit: fsPileup
  
  script:
  prefixfmcool = getPrefix(fmcool)
  prefixfbed = "${fbed.baseName}"
  outPrefix = "${prefixfmcool}.${fbed}.gDist${params.gDistNBins}.ex${excludeZeroSub}"
  options = "--gDistNBins ${params.gDistNBins} --excludeZeroSub ${excludeZeroSub} --cBounds ${params.cvmin} ${params.cvmax} "
  if(params.bClog) {
    options = "${options} --cLog "
  }
  if(bMedian) {
    outPrefix = "${outPrefix}.median"
    options = "${options} --bMedian "
  } else {
    outPrefix = "${outPrefix}.mean"
  }
  if(bNormDiag) {
    outPrefix = "${outPrefix}.normed"
    options = "${options} --bNormDiag "
  }
  if(bIntra) {
    outPrefix = "${outPrefix}.wintra"
    options = "${options} --bIntra "
  }
  if(bInter) {
    options = "${options} --bInter "
  } else {
    outPrefix = "${outPrefix}.nointer"
  }
  outCoo = "${outPrefix}.${params.nBins}x${params.binSize}.mat.txt"
  outSVG = "${outCoo}.svg"
  titleStr = getMicroCSampleName("${fmcool.getSimpleName()}")
  title = "\$\'${titleStr}\'"
  """
  pileupContacts.py ${fmcool} ${fbed} ${params.nBins} ${outCoo} --binSize ${params.binSize} --title ${title} ${options}
  """
}

workflow {
  fmcools = Channel.fromPath(params.inputCools).map { [it.getSimpleName(), it] }

  if(!params.bCmpPeaks2Ways) {
    fPeaks = Channel.fromPath(params.inputPeaks).map { [it.getSimpleName(), it] }
  } else {
    //in the 2-way mode, split the peak file name by the keyword '.up.' to get the two
    //condition names, each of which is one sample type of the micro-c data
    regexSample = /^([^.]+)\.\S+\.([^.]+)\..*$/
    /* fPeaks = Channel.fromPath(params.inputPeaks) */
    Channel.fromPath(params.inputPeaks).multiMap { it ->
        //Get the sample 1 from the first name
        sample1 : [(it.getName() =~ regexSample)[0][1], it]
        //Get the sample 2 from the first name
        sample2 : [(it.getName() =~ regexSample)[0][2], it]
    }.set { fPeakPairs }
    fPeaks = fPeakPairs.sample1.concat(fPeakPairs.sample2)
  }
 
  // match up the sample by name between ChIP peaks and micro-c
  if(params.bMatchCoolsAndPeaks) {
    fmcoolPeakPairs = fmcools.combine(fPeaks, by: 0).map { [it[1], it[2]] }
  } else {
    fmcoolPeakPairs = fmcools.map { it[1] }.combine(fPeaks.map{ it[1] })
  }

  // test different combination of options
  excludeZeroSubs = Channel.from([0])
  bMedians = Channel.from([true])
  bNormDiags = Channel.from([false])
  bIntras = Channel.from([false])
  bInters = Channel.from([true])

  finputs = fmcoolPeakPairs.combine(excludeZeroSubs).combine(bMedians).combine(bNormDiags).combine(bIntras).combine(bInters)
  pileupContacts(finputs)
}
