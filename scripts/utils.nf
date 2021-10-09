#!/usr/bin/env nextflow

process cat {

  cache 'lenient'

  publishDir "${outDir}", mode: 'copy', overwrite: true

  input:
  tuple val(outDir), val(fOut), file(fIn1), file(fIn2)

  output:
  tuple val(outDir), file(fOut), emit:fsCat

  script:
  """
  cat ${fIn1} ${fIn2} > ${fOut}
  """
}

process unpbzip2 {

  executor = 'slurm'
  queue = 'campus-new'
  memory = '22 GB'
  cpus = 8

  cache 'lenient'
 
  input:
  tuple val(outDir), file(fZipped)

  output:
  tuple val("${outDir}"), file(fout), emit: fsUnzipped

  script:
  fout = "${fZipped.getBaseName()}"
  """
  pbzip2 -d ${fZipped} 
  """
}

process pbzip2 {

  executor = 'slurm'
  queue = 'campus-new'
  memory = '22 GB'
  cpus = 8

  cache 'lenient'

  publishDir "${outDir}", mode: 'copy', overwrite: true
 
  input:
  tuple val(outDir), file(f)

  output:
  tuple val(outDir), file(fout), emit: fsBZ2

  script:
  fout = "${f}.bz2"
  """
  pbzip2 ${f} 
  """
}

process sha256 {

  executor = 'slurm'
  queue = 'campus-new'
  memory = '22 GB'
  cpus = 8

  cache 'lenient'

  publishDir "${outDir}", mode: 'copy', overwrite: true

  input:
  tuple outDir, file(f)

  output:
  tuple val(outDir), file(fout), emit: fsSHA256

  script:
  fout = "${f.getName()}.sha256"
  """
  sha256sum ${f} > ${fout}
  """
}

process sha256Check {

  input:
  tuple outDir, file(f), file(fChecksum)

  script:
  """
  sha256sum -c ${fChecksum}
  """
}

params.fchrSizes = "/fh/fast/tsukiyama_t/grp/tsukiyamalab/Sarah/Dejun/reference_genome/sacCer3.chrSizes"
params.mapQ = 6

process pairSam {
  executor = 'slurm'
  queue = 'campus-new'
  memory = '22 GB'
  cpus = 8

  cache 'lenient'
  publishDir "${outDir}", mode: 'copy', overwrite: true

  input:
  tuple outDir, file(fSam1), file(fSam2), nINCutoff 

  output:
  tuple outDir, file(fPairIN), file(fPairOUT), file(fPairSAME), file(fPair), emit: fsPair 
  file flog

  script:
  regexp = /^(.*).R[12].*$/
  fPrefix = (fSam1.getName() =~ regexp)[0][1]
  fPrefixPair = "${fPrefix}"
  if( nINCutoff == 0 ) {
    fPrefixPair += ".NoFilter"
  } else if( nINCutoff == 9999999999 ) {
    fPrefixPair += ".NoIN"
  } else {
    fPrefixPair += ".FilterIN${nINCutoff}"
  }
  //This is the log file from filter_and_sort_pairs
  flog = "${fPrefixPair}.filtersort.log"
  fPairIN = "${fPrefixPair}.IN.pair"
  fPairOUT = "${fPrefixPair}.OUT.pair"
  fPairSAME = "${fPrefixPair}.SAME.pair"
  fPair = "${fPrefixPair}.pair"

  //These are filtered according to input cutoff of IN reads
  //The first call should generate 3 files and the 
  //other 2 calls verify the sha256
  //Note that Nextflow only allows 1 shell command block so 
  //the commands from different blocks have to be concatenated
  """
  filter_and_sort_pairs ${fSam1} ${fSam2} ${params.fchrSizes} ${params.mapQ} ${nINCutoff} 1 ${fPrefixPair} && \
  filter_and_sort_pairs ${fSam1} ${fSam2} ${params.fchrSizes} ${params.mapQ} ${nINCutoff} 0 ${fPrefixPair}
  """
}

process gzip {
  executor = 'slurm'
  queue = 'campus-new'
  memory = '22 GB'
  cpus = 8

  cache 'lenient'
  publishDir "${outDir}", mode: 'copy', overwrite: true

  input:
  tuple outDir, file(f)

  output:
  tuple outDir, file(fout), emit: fsZipped

  script:
  fout = "${f.getName()}.gz"
  """
  pbgzip ${f}
  """
}

params.juicerJar = "~/local/juicer_tools_1.12.03.jar"

process pair2hic {
  executor = 'slurm'
  queue = 'campus-new'
  memory = '32 GB'
  cpus = 8

  cache 'lenient'
  publishDir "${outDir}", mode: 'copy', overwrite: true

  input:
  tuple outDir, file(fpairGZ)

  output:
  tuple outDir, file(fout), emit: fsHic

  script:
  fout = "${fpairGZ.getName()}.hic"
  """
  java -Xmx32g -jar ${params.juicerJar} pre -r 10,50,100,150,200,500,1000,3200,5000 ${fpairGZ} ${fout} sacCer3
  """
}

process hic2mcool {
  executor = 'slurm'
  queue = 'campus-new'
  memory = '22 GB'
  cpus = 8

  cache 'lenient'
  publishDir "${outDir}", mode: 'copy', overwrite: true

  input:
  tuple outDir, file(fhic)

  output:
  tuple outDir, file(fout), emit: fsMcool

  script:
  fout = "${fhic.getName()}.mcool"
  """
  if command -v pyenv > /dev/null; then eval "\$(pyenv init -)"; fi && \
  pyenv shell hic2cool && \
  hic2cool convert ${fhic} ${fout}
  """
}

def getMicroCSampleName ( samplePrefix ) {
  def nameMap = [
    "Q_from_paper" : "Q",
    "L_from_paper" : "Log",
    "Q_5toA_merged" : "Q 5toA",
    "Q_R17R19A_merged" : "Q R17R19A",
    "Q_TSA_merged" : "Q TSA",
    "Q_dCondensin_6909" : "Q -condensin",
    "Q_XS_merged" : "Q XS",
    "L_XS_merged" : "Log XS",
    "Q_dNhp6A_7334" : "Q -Nhp6A",
    "Q_HHF2_7206" : "Q HHF2",
  ]
  return nameMap[samplePrefix]
}


def getPrefix( file ) {
  regexpPE = /^(.*).pair.*$/
  //Per https://www.nextflow.io/docs/latest/script.html#capturing-groups
  //match[0] is [whole-string, group1, group2, ...]
  //match[0][1] is group1
  (file.baseName =~ regexpPE)[0][1]
}

def splitPrefix ( file ) {
  re = /^(.*)-(.*)[.]\d+.*$/
  m = (file.getName() =~ re)
  [m[0][1], m[0][2]]
}

def trimSuffix ( s ) {
  s.take(s.lastIndexOf('.'))
}

def getSampleName ( file ) {
  // get sample name assuming the naming scheme:
  // [celltype]_[treatment]_[replicaIndex]
  // so that [celltype]_[treatment] is returned
  // as id of sample name
  s = file.simpleName
  return s.take(s.lastIndexOf('_'))
}

process contactDecay {
  executor = 'slurm'
  queue = 'campus-new'
  memory = '22 GB'
  cpus = 8
  cache = 'lenient'
  errorStrategy { 
    task.exitStatus == 9 ||
    task.exitStatus == 64 ||
    task.exitStatus == 125 ||
    task.exitStatus == 130 ||
    task.exitStatus == 131 ||
    task.exitStatus == 134 ||
    task.exitStatus == 137 ||
    task.exitStatus == 139 ||
    task.exitStatus == 140 ? 
    'terminate' : 'retry' }
  maxRetries = 9999
  publishDir "${outdir}", mode: 'copy', overwrite: true

  input:
  tuple file(f), file(fbed), val(outdir), val(binSize), val(bNormNpairs), val(bSgrNucPeaks), val(bGroupByGeneCluster) 

  output:
  file "${foutprefix}*"

  script:
  foutprefix = "${getPrefix(f)}.${binSize}bp"
  options = ""
  if ( bSgrNucPeaks ) {
    sampleName = f.getSimpleName()
    options = "${options} --fNucSgr /fh/fast/tsukiyama_t/grp/tsukiyamalab/Sarah/Dejun/nucleosome_peaks/${sampleName}.sgr"
    foutprefix += ".sgr"
  }
  if ( bNormNpairs ) {
    options = "${options} --bNormNpairs"
    foutprefix += ".normedNpairs"
  }
  if ( bGroupByGeneCluster ) {
    //The 12th column of a bed file is the cluster index
    options = "${options} --fbed ${fbed} --col 12" 
  }
  //foutWholeGenome= "${foutprefix}.WholeGenome.contactDecay.mean"
  """
  #!/bin/bash
  set -euxo pipefail
  contactDecay.py ${f} ${binSize} ${foutprefix} ${options}
  """
}

process mergePairGZ {
  executor = 'slurm'
  queue = 'campus-new'
  memory = '22 GB'
  cpus = 8
  cache = 'lenient'
  errorStrategy { 
    task.exitStatus == 9 ||
    task.exitStatus == 64 ||
    task.exitStatus == 125 ||
    task.exitStatus == 130 ||
    task.exitStatus == 131 ||
    task.exitStatus == 134 ||
    task.exitStatus == 137 ||
    task.exitStatus == 139 ||
    task.exitStatus == 140 ? 
    'terminate' : 'retry' }
  maxRetries = 9999
  publishDir "${outDir}", mode: 'copy', overwrite: true

  input:
  tuple outDir, fOutPrefix, file(f1), file(f2)

  output:
  tuple val(outDir), file(fout), emit: fsMerged

  script:
  fout = "${fOutPrefix}.pair.gz"
  //Needs to clone https://github.com/4dn-dcic/pairix and install it
  """
  ~/Downloads/pairix/util/merge-pairs.sh ${fOutPrefix} ${f1} ${f2} && \
  mv ${fOutPrefix}.pairs.gz ${fout}
  """
}

def file2path ( f ) {
  path = file(f).toRealPath()
  dir = path.getParent()
  return [ dir, file(f) ]
}

params.genomeIndex = "/fh/fast/tsukiyama_t/grp/tsukiyamalab/Sarah/Dejun/reference_genome/sacCer3/bowtie2/sacCer3"

process bowtie2 {

  executor = 'slurm'
  queue = 'campus-new'
  memory = '44 GB'
  cpus = 16

  cache 'lenient'
 
  input:
  tuple val(outDir), file(fFastq)

  output:
  tuple val(outDir), file(fout), emit: fsSAM

  script:
  fout = "${fFastq}.bowtie2.sam"
  flog = "${fFastq}.bowtie2.log"
  """
  bowtie2 --very-sensitive -t -p ${task.cpus} --reorder -x ${params.genomeIndex} -U ${fFastq} -S ${fout} 2> ${flog}  
  """
}


params.genes = "/fh/fast/tsukiyama_t/grp/tsukiyamalab/Sarah/Dejun/Sc.cerevisiae.annotation_Steinmetz_2013.gtf"


process computeMatrix {

  executor = 'slurm'
  queue = 'campus-new'
  memory = '22 GB'
  cpus = 8

  cache 'lenient'

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  tuple file(fBigwig), val(mode), val(a), val(b)

  output:
  file(fMatrixGZ)

  script:
  fMatrixGZ = "${fBigwig}.TSS"
  if(mode == "reference-point") {
    fMatrixGZ = "${fMatrixGZ}a${a}b${b}"
  } else if(mode == "scale-regions") {
    fMatrixGZ = "${fMatrixGZ}TES"
  }
  fMatrixGZ = "${fMatrixGZ}.matrix.gz"
  """
  computeMatrix ${mode} -a ${a} -b ${b} -R ${params.genes} -S ${fBigwig} -o ${fMatrixGZ}
  """
 
}

process plotHeatmapDeepTools {

  executor = 'slurm'
  queue = 'campus-new'
  memory = '22 GB'
  cpus = 8

  cache 'lenient'

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  tuple file(fMatrixGZ), val(vMin), val(vMax)

  output:
  file(fSVG)

  script:
  fSVG = "${fMatrixGZ.getName()}.svg"
  """
  plotHeatmap -m ${fMatrixGZ} -out ${fSVG} --regionsLabel "Genes" --plotTitle ${fMatrixGZ}  --colorMap 'plasma' --sortUsing sum --zMin ${vMin} --zMax ${vMax} --yMin ${vMin} --yMax ${vMax}  --averageTypeSummaryPlot median
  """
}

process diffDenseMatrix {
  executor = 'slurm'
  queue = 'campus-new'
  memory = '22 GB'
  cpus = 8

  cache 'lenient'

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  tuple file(f1), file(f2), val(fout)

  output:
  file(fout)

  script:
  """
  matrixDiff.py ${f1} ${f2} ${fout}
  """
}

process plotHeatmap {
  executor = 'slurm'
  queue = 'campus-new'
  memory = '22 GB'
  cpus = 8

  cache 'lenient'

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  tuple file(fin), val(fout), val(title)

  output:
  file(fout)
  
  script:
  options = "--bTicks0OriginCenter --bTicks1OriginCenter --tickLabelMultiplier ${params.tickLabelMultiplier} --tickEvery ${params.tickEvery} --title ${title} --cLabel ${params.cLabel} --cBounds ${params.cvmin} ${params.cvmax} --cMap ${params.cmap} --xLabel ${params.xLabel} --yLabel ${params.yLabel} "
  if(params.bClog) {
    options = "${options} --cLog "
  }
  """
  plotHeatmap.py ${fin} ${fout} ${options}
  """
}

