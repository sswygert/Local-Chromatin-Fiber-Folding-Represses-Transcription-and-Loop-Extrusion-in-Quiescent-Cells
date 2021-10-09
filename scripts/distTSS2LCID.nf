#!/usr/bin/env nextflow
nextflow.preview.dsl=2

params.inputTSS = "./transcripts/*.txt"
params.inputLCID = "./LCID/Q_from_paper.Brn10.LCID.MACSAll.bed"
params.outdir = "./distTSS2LCID"

include { getPrefix } from './utils'

process distTSS2LCID {
  executor = 'slurm'
  queue = 'campus-new'
  memory = '22 GB'
  cpus = 8

  cache 'lenient'

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  tuple file(fTSS), file(fLCID)

  output:
  file(fout)

  script:
  fout = "${fTSS.getBaseName()}.dist.${fLCID.getBaseName()}.txt"
  """
  distTSS2LCID.py ${fTSS} ${fLCID} ${fout}
  """
}

workflow {
  fTSSs = Channel.fromPath(params.inputTSS)
  fLCIDs = Channel.fromPath(params.inputLCID)
  fInputs = fTSSs.combine(fLCIDs)
  distTSS2LCID(fInputs)
}
