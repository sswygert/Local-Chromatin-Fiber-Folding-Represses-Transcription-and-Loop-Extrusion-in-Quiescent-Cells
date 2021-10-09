#!/usr/bin/env nextflow
nextflow.preview.dsl=2
/* include './utils' */
include { file2path; sha256; hic2mcool } from './utils'

params.input = "./micro-c/Q_XS_5783Replicate/*.hic"

workflow {
  fsHic = Channel.fromPath(params.input).map { file2path(it) }
  hic2mcool(fsHic)
  sha256(hic2mcool.out.fsMcool)
}
