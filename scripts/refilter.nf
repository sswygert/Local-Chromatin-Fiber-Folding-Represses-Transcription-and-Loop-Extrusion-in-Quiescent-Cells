#!/usr/bin/env nextflow
nextflow.preview.dsl=2
/* include './utils' */
include { unpbzip2; sha256; gzip; pairSam; pairSam as pairSamFilterIN; pair2hic; hic2mcool } from './utils'

params.input = "./micro-c/*/*.sam.bz2"

def file2path ( f ) {
  path = file(f).toRealPath()
  dir = path.getParent()
  return [ dir, file(f) ]
}

workflow {
  fsSamZipped = Channel.fromPath(params.input).map { file2path(it) }
  unpbzip2(fsSamZipped)
  fsSam = unpbzip2.out.fsUnzipped
  fsSamPair = fsSam.groupTuple().map { [it[0], 
    //make sure R1 goes before R2
    it[1][0].getName() < it[1][1].getName() ? it[1][0] : it[1][1],
    it[1][0].getName() < it[1][1].getName() ? it[1][1] : it[1][0],
    ] }
  //No filtering
  fsSamPairNoFilter = fsSamPair.combine(Channel.from(0))
  fsSamPairFilterIN = fsSamPair.combine(Channel.from(9999999999))
  pairSam(fsSamPairNoFilter)
  pairSamFilterIN(fsSamPairFilterIN)
  fsPairNoFilter = pairSam.out.fsPair
  fsPairFilterIN = pairSamFilterIN.out.fsPair
  fsPairAll = fsPairNoFilter.mix(fsPairFilterIN).
  //group files from same dir
  groupTuple().
  //one tuple for each dir ->  each file has its own associated dir
  flatMap {
    ftuple -> [ 
      //outDir   //fPair
      ftuple[0], ftuple[1][0], //NoFilter IN
      ftuple[0], ftuple[1][1], //NoIn IN
      ftuple[0], ftuple[2][0], //... OUT
      ftuple[0], ftuple[2][1], //... OUT
      ftuple[0], ftuple[3][0], //... SAME
      ftuple[0], ftuple[3][1], //... SAME
      ftuple[0], ftuple[4][0], //... Merged
      ftuple[0], ftuple[4][1], //... Merged
      ]
  }.collate(2)  gzip(fsPairAll)
  pair2hic(gzip.out.fsZipped)
  hic2mcool(pair2hic.out.fsHic)
  sha256(hic2mcool.out.fsMcool)
}
