#!/usr/bin/env nextflow
nextflow.preview.dsl=2
include { file2path; cat; unpbzip2; sha256; pair2hic; hic2mcool; mergePairGZ } from './utils'
include { fastq2mcool } from './workflows'

def getSuffix( file ) {
  regexpPE = /^[^.]+([.].*)[.]pair[.]gz$/
  //Per https://www.nextflow.io/docs/latest/script.html#capturing-groups
  //match[0] is [whole-string, group1, group2, ...]
  //match[0][1] is group1
  (file =~ regexpPE)[0][1]
}

def getPrefix( file ) {
  regexpPE = /^(.*).pair.*$/
  //Per https://www.nextflow.io/docs/latest/script.html#capturing-groups
  //match[0] is [whole-string, group1, group2, ...]
  //match[0][1] is group1
  (file.baseName =~ regexpPE)[0][1]
}

def getReadIdx( file ) {
  regexPE = /(R[12]).fastq/
  return (file =~ regexPE)[0][1]
}

//~/Downloads/pairix/util/merge-pairs.sh  L_XS_merged.NoIN.IN ../L_XS_5781/L_XS_5781.NoIN.IN.pair.gz ../L_XS_5783/L_XS_5783.NoIN.IN.pair.gz

params.dataDir = "./micro-c"
params.Rep1 = ""
params.Rep2 = ""
params.out = ""
params.bFASTQ = true //if merge fastq and start from scratch

outDir = "${params.dataDir}/${params.out}"

workflow {

  if(params.bFASTQ) {
    input1 = "${params.dataDir}/${params.Rep1}/*.fastq.bz2"
    input2 = "${params.dataDir}/${params.Rep2}/*.fastq.bz2"
    fs1FASTQBZ2 = Channel.fromPath(input1)
    fs2FASTQBZ2 = Channel.fromPath(input2)

    //Here I assume there are only two: R1 and R2 fastq files for each replica
    fsFASTQBZ2 = fs1FASTQBZ2.mix(fs2FASTQBZ2).map {
      [getReadIdx(it), it]
    }.groupTuple().map {
      [
        outDir,
        "${params.out}.${it[0]}.fastq.bz2",
        it[1][0].getName() ==~ /${params.Rep1}/ ? it[1][0] : it[1][1],
        it[1][0].getName() ==~ /${params.Rep1}/ ? it[1][1] : it[1][0],
      ]
    }
    cat(fsFASTQBZ2)
    unpbzip2(cat.out.fsCat)
    fastq2mcool(unpbzip2.out.fsUnzipped)

  } else {
    input1 = "${params.dataDir}/${params.Rep1}/*.pair.gz"
    input2 = "${params.dataDir}/${params.Rep2}/*.pair.gz"
    fsPairRep1 = Channel.fromPath(input1)
    fsPairRep2 = Channel.fromPath(input2)
    
    fsPairs = fsPairRep1.mix(fsPairRep2).map {
      [getSuffix(it), it]
    }.groupTuple().map { 
      [
        outDir,
        "${params.out}${it[0]}",
        it[1][0].getName() < it[1][1].getName() ? it[1][0] : it[1][1],
        it[1][0].getName() < it[1][1].getName() ? it[1][1] : it[1][0],
      ]
    }
    
    mergePairGZ(fsPairs)
    fsPairMerged = mergePairGZ.out.fsMerged
    pair2hic(fsPairMerged)
    hic2mcool(pair2hic.out.fsHic)
    sha256(hic2mcool.out.fsMcool)
  }

}
