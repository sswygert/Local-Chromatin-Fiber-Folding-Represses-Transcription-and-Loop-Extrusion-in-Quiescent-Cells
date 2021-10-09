#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include { getMicroCSampleName } from './utils'

params.input = "./contactDecay/{Q_from_paper,L_from_paper,Q_XS_merged,L_XS_merged,Q_R17R19A_merged,Q_5toA_merged,Q_TSA_merged,Q_dCondensin_6909}*.NoIN.SAME*.10bp*WholeGenome*.mean"
params.outdir = "./contactDecay"
params.min1 = 0
params.max1 = 500
params.min2 = 500
params.max2 = 1000

process contactOddRatio {
  executor = 'slurm'
  queue = 'campus-new'
  memory = '2 GB'
  cpus = 1

  cache 'lenient'

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  tuple file(fin), val(min1), val(max1), val(min2), val(max2)

  output:
  tuple file(fout), stdout, emit: oddRatio
  
  script:
  fout = "${fin}.oddRatio.${min1}-${max1}.vs.${min2}-${max2}.txt"
  """
  contactOddRatio.py ${fin} ${fout} ${min1} ${max1} ${min2} ${max2}
  tail -n1 ${fout}

  """
}

workflow {
  fsinput = Channel.fromPath(params.input)
  input = fsinput.map {[it, params.min1, params.max1, params.min2, params.max2]} 

  contactOddRatio(input)

  contactOddRatio.out.oddRatio.map{
    "${getMicroCSampleName(it[0].getSimpleName()).replaceAll("\\s+", "_")} ${it[1]}"}.
    collectFile(
      storeDir: "${params.outdir}",
      name: "contactOddRatio.txt",
      seed: "#Sample p(${params.min1} <= d <= ${params.max1}) p(${params.min2} <= d <= ${params.max2}) p2/p1\n"
      )
}
