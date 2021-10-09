#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include { file2path; pbzip2 } from './utils'

params.input = "./micro-c/*/*.fastq"
params.outdir = "./micro-c"

workflow {
  fs = Channel.fromPath(params.input).map { file2path(it) }
  pbzip2(fs)
}
