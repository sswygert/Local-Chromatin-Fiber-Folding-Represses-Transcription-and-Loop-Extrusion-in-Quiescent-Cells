#!/usr/bin/env nextflow
nextflow.preview.dsl=2

params.input = "./micro-c/{Q_HHF2_7177ReplicateResequenced,Q_TSA_5781Resequenced}/*.fastq.gz"
params.outdir = "./micro-c"

process zcatFastqGZ {

  executor = 'slurm'
  queue = 'campus-new'
  memory = '22 GB'
  cpus = 8

  cache 'lenient'

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  tuple val(foutPrefix), file(fsgz)

  output:
  file(fout)

  script:
  fout = "${foutPrefix}.fastq"
  """
  zcat ${fsgz} > ${fout}
  """
}

workflow {

  fs = Channel.fromPath(params.input)
       .map { file ->
         def path = file.getParent()
         def prefix = path.getName()
         readTagRegex = /_(R[12])_/
         def readTag = (file.name.toString() =~ readTagRegex)[0][1]
         def foutPrefix = "${prefix.toString()}_${readTag}"
         return tuple(foutPrefix, file)
       }
       .groupTuple()
       .map { it ->
         def foutPrefix = it[0]
         def fs = it[1]
         fsSorted = fs.sort { (it =~ /_([0-9]+).fastq.gz$/)[0][1]  }
         return tuple(foutPrefix, fsSorted)
       }

  zcatFastqGZ(fs)

}


