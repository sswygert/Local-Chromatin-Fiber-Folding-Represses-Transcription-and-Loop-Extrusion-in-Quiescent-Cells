#!/usr/bin/env nextflow
nextflow.preview.dsl=2

params.input = "./LCID_micro-c_pileup/*All*median*nointer*.txt"
params.outdir = "./LCID_micro-c_pileup"
params.cvmin = -100
params.cvmax = 100
params.bClog = true
params.cLabel = "\'Difference in median counts\'"
params.cmap = "bwr"
// this sets the x- and y-ticks in kb assuming 200 bp bin size in the heatmap
params.tickLabelMultiplier = 200 / 1000
// this plot every 2 kb assuming 200 bp bin size
params.tickEvery = 2000 / 200
params.xLabel = "\'Genomic distance (kb)\'"
params.yLabel = "\'Genomic distance (kb)\'"

include { diffDenseMatrix; plotHeatmap; getMicroCSampleName } from './utils'

workflow {
  fmats = Channel.fromPath(params.input)
  fmatPairs = fmats.combine(fmats).filter { it[0] != it[1] }
  fDiffIns = fmatPairs.map { [it[0], it[1], "${it[1].getBaseName()}-${it[0].getBaseName()}.txt"] }
  diffDenseMatrix(fDiffIns)
  fDiffMats = diffDenseMatrix.out.map { [ "${it.getBaseName()}", it ] }
  fSVGs = fmatPairs.map { ["${it[1].getBaseName()}-${it[0].getBaseName()}", "${it[1].getBaseName()}-${it[0].getBaseName()}.svg"] }
  fTitles = fmatPairs.map { ["${it[1].getBaseName()}-${it[0].getBaseName()}", "\$\'${getMicroCSampleName(it[1].getSimpleName())} - ${getMicroCSampleName(it[0].getSimpleName())}\'"] }
  fsInHeatmap = fDiffMats.join(fSVGs, by: 0).join(fTitles, by: 0).map { [it[1], it[2], it[3]] }
  plotHeatmap(fsInHeatmap)
}
