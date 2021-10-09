#!/usr/bin/env nextflow
nextflow.preview.dsl=2
/* include './utils' */
include { file2path; pbzip2 } from './utils'
include { fastq2mcool } from './workflows'

// params.input = "./micro-c/*/*.fastq"
params.input = "./micro-c_archive/Q_HHF2_7177ReplicateResequenced/*.fastq"

workflow {
  fsFastq = Channel.fromPath(params.input).map { file2path(it) }
  fastq2mcool(fsFastq)
  //pbzip2 the fastq and sam files
  pbzip2(fsFastq)
}
