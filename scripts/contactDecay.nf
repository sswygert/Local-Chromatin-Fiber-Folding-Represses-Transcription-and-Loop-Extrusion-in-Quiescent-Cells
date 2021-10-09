#!/usr/bin/env nextflow
nextflow.preview.dsl=2
include { contactDecay } from './utils'

//examine the contact count ratio between IN, OUT and SAME orientations
params.outdir = "./contactDecay"
//params.fsmcool = "./micro-c/*/*.NoIN*.pair*.mcool"
params.fsmcool = "./micro-c/{Q_from_paper,L_from_paper}/*.NoIN.SAME*.pair*.mcool"
params.fsbed = "./gene_clusters/{Q_from_paper,L_from_paper}*2clust.bed"
params.binSize = 10
params.bNormNpairs = true
params.bSgrNucPeaks = false
params.bGroupByGeneCluster = true

workflow {
  fsmcool = Channel.fromPath(params.fsmcool).map { [it.getSimpleName(), it] }
  fsbed = Channel.fromPath(params.fsbed).map { [it.getSimpleName(), it] }

  //match mcool and bed file by sample name
  fsinput = fsmcool.combine(fsbed, by: 0).map{ [it[1], it[2]] }.\
    combine(Channel.from(params.outdir)).\
    combine(Channel.from(params.binSize)).\
    combine(Channel.from(params.bNormNpairs)).\
    combine(Channel.from(params.bSgrNucPeaks)).\
    combine(Channel.from(params.bGroupByGeneCluster))
  
  contactDecay(fsinput)
  fsCD = contactDecay.out
}
