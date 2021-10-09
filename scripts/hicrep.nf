#!/usr/bin/env nextflow
nextflow.preview.dsl=2


params.tag = "NoIN.pair"
/* // for intra-replica hicrep */
/* params.input = "./micro-c_archive/{L_XS_5781,L_XS_5783,Q_R17R19A_7200,Q_R17R19A_7207,Q_5toA_7175,Q_5toA_7208,Q_TSA_5783Merged,Q_TSA_5781Merged,Q_XS_5783Merged,Q_XS_5781}/*.${params.tag}.*.mcool" */
/* params.bPerSample = true */
/* params.outdir = "./hicrep/perSample" */

// for inter merged hicrep
params.input = "./micro-c/{L_from_paper,L_XS_merged,Q_from_paper,Q_5toA_merged,Q_R17R19A_merged,Q_HHF2_7206,Q_TSA_merged,Q_XS_merged,Q_dCondensin_6909}/*.${params.tag}.*.mcool"
params.bPerSample = false
params.outdir = "./hicrep/mergedSample"

//whether to keep the pattern {prefix}.{tag} as in {prefix}.{tag}.pair.* as
//output pattern instead of {prefix} only 
params.bKeepSuffix = false
params.binSize = 5000

include { getPrefix; splitPrefix; trimSuffix; getSampleName } from './utils'

process HiCRep {
  executor = 'slurm'
  queue = 'campus-new'
  memory = '22 GB'
  cpus = 8

  cache 'lenient'

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  tuple file(f1), file(f2)

  output:
  path(foutMedian, emit: fSCCSMedian)
  path(fout, emit: fSCCS)

  script:
  p1=getPrefix(f1)
  p2=getPrefix(f2)
  fout="${p1}-${p2}.${params.binSize}.sccs.txt"
  foutMedian="${p1}-${p2}.${params.binSize}.sccs.txt.median"
  n1=p1
  n2=p2
  if (!params.bKeepSuffix) {
    n1 = trimSuffix(n1)
    n2 = trimSuffix(n2)
  }

  shell:
  '''
  set -euxo pipefail
  hicrep.py !{f1} !{f2} !{fout} --binSize !{params.binSize} --h 1 --dBPMax 100000 && \
  median=$(python3 -c "import numpy as np; print('%30.15e' % np.median(np.loadtxt('!{fout}').astype(float)))") && \
  echo !{n1} !{n2} ${median} > !{foutMedian} 

  '''
}

process selfHiCRep {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  file(f)

  output:
  file(fout)

  script:
  p = getPrefix(f)
  fout="${p}-${p}.${params.binSize}.sccs.txt.median"
  n = p
  if (!params.bKeepSuffix) {
    n = trimSuffix(n)
  }

  shell:
  '''
  echo !{n} !{n} 1.0 > !{fout} 
  '''
}

workflow {
  fs = Channel.fromPath(params.input)
  fPairs = fs.combine(fs).filter { it[0] < it[1] }
  if(params.bPerSample) {
    // only do intra-sample hicrep
    fPairs = fPairs.filter { getSampleName(it[0]) == getSampleName(it[1]) }
  }
  HiCRep(fPairs)
  selfHiCRep(fs)
  fAllSCCS = HiCRep.out.fSCCSMedian.concat(selfHiCRep.out)
  fmerged = "${params.outdir}/${params.tag}.${params.binSize}.txt"
  fAllSCCS.collectFile(name: fmerged, seed:"sample1 sample2 scc\n")
}
