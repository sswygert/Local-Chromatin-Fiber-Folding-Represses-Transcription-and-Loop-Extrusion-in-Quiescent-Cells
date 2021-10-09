#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/* params.input = "./micro-c/{Q_from_paper,L_from_paper,Q_5toA_merged,Q_R17R19A_merged,Q_dCondensin_6909,Q_dNhp6A_7334,Q_XS_merged,L_XS_merged}/*.NoIN.{OUT,SAME}.*.mcool" */
params.input = "./micro-c/{Q_from_paper,L_from_paper}/*.NoIN.{OUT,SAME}.*.mcool"
params.outdir = "./np2np1ratio"
params.bNucPeak = true

/* include './utils' */
include { computeMatrix; computeMatrix as computeMatrixTSSTES; 
  computeMatrix as computeMatrixRatioDiffSum; computeMatrix as computeMatrixRatioDiffSumTSSTES; 
  plotHeatmapDeepTools; plotHeatmapDeepTools as plotHeatmapDeepToolsTSSTES;
  plotHeatmapDeepTools as plotHeatmapDeepToolsRatioDiffSum; plotHeatmapDeepTools as plotHeatmapDeepToolsRatioDiffSumTSSTES;
  } from './utils'
include { getPrefix } from './utils'


process np2np1ratio {
  executor = 'slurm'
  queue = 'campus-new'
  memory = '22 GB'
  cpus = 8

  cache 'lenient'

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  file(fmcool) 

  output:
  tuple file(outNp1), file(outNp2), file(outRatio), file(outRatioNormed), file(outRatioDiffSum), emit: fsRatio
  
  script:
  outPrefix = getPrefix(fmcool)

  outNp1 = "${outPrefix}.np1.bw"
  outNp2 = "${outPrefix}.np2.bw"
  outRatio = "${outPrefix}.np2np1.bw"
  outRatioNormed = "${outPrefix}.np2np1normed.bw"
  outRatioDiffSum = "${outPrefix}.np2np1diffsum.bw"


  options = "${fmcool} 50 1000 "

  if(params.bNucPeak) {
    //use the nucleosome peak file
    sampleName = fmcool.getSimpleName()
    //use the nucleosome peak files to run np2np1ratio.py
    outPrefix = "${outPrefix}.sgr"
    outNp1 = "${outPrefix}.np1.bw"
    outNp2 = "${outPrefix}.np2.bw"
    outRatio = "${outPrefix}.np2np1.bw"
    outRatioNormed = "${outPrefix}.np2np1normed.bw"
    outRatioDiffSum = "${outPrefix}.np2np1diffsum.bw"
    options = "${options} ${outPrefix} --fNucSgr /fh/fast/tsukiyama_t/grp/tsukiyamalab/Sarah/Dejun/nucleosome_peaks/${sampleName}.sgr"
  } else {
    //Parameters for SAME orientations
    np1bpBegin = 50
    np1bpLast = 250
    np2bpBegin = 250
    np2bpLast = 450
    if(outPrefix ==~ /[.]OUT$/) {
      //Parameters for OUT orientations
      np1bpBegin = 0
      np1bpLast = 200
      np2bpBegin = 200
      np2bpLast = 400
    }
    options = "${options} ${outPrefix} --np2np1Bins ${np1bpBegin} ${np1bpLast} ${np2bpBegin} ${np2bpLast}"
  }

  """
  np2np1ratio.py ${options}
  """

}

workflow {
  fs = Channel.fromPath(params.input)
  //compute the ratio
  np2np1ratio(fs)

  //get the np2np1 raw ratio
  fsRatioBigWig = np2np1ratio.out.fsRatio.map { it[2] }
  //generate gene-by-ratio matrix
  //For TSS +/- 3K
  fsRatioBigWigTSS3K = fsRatioBigWig.map { [it,  "reference-point", 3000, 3000] }
  computeMatrix(fsRatioBigWigTSS3K)
  fsMatrixGZTSS3K = computeMatrix.out 
  //For TSS-TES
  fsRatioBigWigTSSTES = fsRatioBigWig.map { [it,  "scale-regions", 0, 0] }
  computeMatrixTSSTES(fsRatioBigWigTSSTES)
  fsMatrixGZTSSTES = computeMatrixTSSTES.out
  //Plot matrix heatmap
  plotHeatmapDeepTools(fsMatrixGZTSS3K.map{[it, 0.4, 1.8]})
  plotHeatmapDeepToolsTSSTES(fsMatrixGZTSSTES.map{[it, 0.4, 1.8]})

  //get the (np2-np1) / (np2+np1) ratio
  fsRatioDiffSumBigWig = np2np1ratio.out.fsRatio.map { it[4] }
  //generate gene-by-ratio matrix
  //For TSS +/- 3K
  fsRatioDiffSumBigWigTSS3K = fsRatioDiffSumBigWig.map { [it,  "reference-point", 3000, 3000] }
  computeMatrixRatioDiffSum(fsRatioDiffSumBigWigTSS3K)
  fsMatrixRatioDiffSumGZTSS3K = computeMatrixRatioDiffSum.out 
  //For TSS-TES
  fsRatioDiffSumBigWigTSSTES = fsRatioDiffSumBigWig.map { [it,  "scale-regions", 0, 0] }
  computeMatrixRatioDiffSumTSSTES(fsRatioDiffSumBigWigTSSTES)
  fsMatrixRatioDiffSumGZTSSTES = computeMatrixRatioDiffSumTSSTES.out
  //Plot matrix heatmap
  plotHeatmapDeepToolsRatioDiffSum(fsMatrixRatioDiffSumGZTSS3K.map{[it, -1.0, 1.0]})
  plotHeatmapDeepToolsRatioDiffSumTSSTES(fsMatrixRatioDiffSumGZTSSTES.map{[it, -1.0, 1.0]})
}
