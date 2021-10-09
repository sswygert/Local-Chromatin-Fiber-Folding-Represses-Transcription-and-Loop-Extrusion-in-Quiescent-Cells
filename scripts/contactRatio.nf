#!/usr/bin/env nextflow

//examine the contact count ratio between IN, OUT and SAME orientations
params.indir = "./contactDecay"
params.outdir = "./contactRatio"
fsOUTIN = Channel.fromFilePairs("${params.indir}/*.{OUT,IN}*", flat: true)
fsSAMEIN = Channel.fromFilePairs("${params.indir}/*.{SAME,IN}*", flat: true)
fsOUTSAME = Channel.fromFilePairs("${params.indir}/*.{OUT,SAME}*", flat: true)

process contactRatio {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true, pattern: "*.txt"
  publishDir "${params.outdir}/plots", mode: 'copy', overwrite: true, pattern: "*.svg"

  input:
  set id, fINp1, fOUTp1 from fsOUTIN
  set id, fINp2, fSAMEp2 from fsSAMEIN
  set id, fOUTp3, fSAMEp3 from fsOUTSAME

  output:
  file foutp1 into fsoutp1
  file foutp2 into fsoutp2
  file foutp3 into fsoutp3
  file fplot into fplots

  script:
  foutp1 = "${fOUTp1.getBaseName()}-${fINp1.getBaseName()}.ratio.txt"
  foutp2 = "${fSAMEp2.getBaseName()}-${fINp2.getBaseName()}.ratio.txt"
  foutp3 = "${fOUTp3.getBaseName()}-${fSAMEp3.getBaseName()}.ratio.txt"
  foutprefix = "${fOUTp1.getSimpleName()}"
  fplot = "${foutprefix}.svg"

  shell:
  '''
  #!/usr/bin/env python3
  import numpy as np
  import matplotlib.pyplot as plt
  import matplotlib
  matplotlib.use('Agg')

  def getRatio(f1, f2):
      c1 = np.loadtxt(f1, usecols=(1,), dtype=float)
      c2 = np.loadtxt(f2, usecols=(1,), dtype=float)
      s = min(c1.size, c2.size)
      gdist1 = np.loadtxt(f1, usecols=(0,), dtype=float)[:s]
      gdist2 = np.loadtxt(f2, usecols=(0,), dtype=float)[:s]
      assert (gdist1 == gdist2).all(), "Genomic distance in the two inputs are different"
      # only take the minimal size
      return gdist1, c1[:s] / c2[:s]

  gdist1, ratio1 = getRatio("!{fOUTp1}", "!{fINp1}") 
  gdist2, ratio2 = getRatio("!{fSAMEp2}", "!{fINp2}") 
  gdist3, ratio3 = getRatio("!{fOUTp3}", "!{fSAMEp3}") 
  np.savetxt("!{foutp1}", np.vstack((gdist1,ratio1)).T, "%30.15e", header="contactRatio.nf")
  np.savetxt("!{foutp2}", np.vstack((gdist2,ratio2)).T, "%30.15e", header="contactRatio.nf")
  np.savetxt("!{foutp3}", np.vstack((gdist3,ratio3)).T, "%30.15e", header="contactRatio.nf")
  plt.title("!{foutprefix}")
  plt.xlim(0, 2000)
  plt.yscale('log', basey=2)
  plt.ylim(2**(-10), 2**9)
  plt.plot(gdist1, ratio1, label="OUT/IN")
  plt.plot(gdist2, ratio2, label="SAME/IN")
  plt.plot(gdist3, ratio3, label="OUT/SAME")
  plt.legend()
  plt.xlabel("Genomic Distance (10 bp)")
  plt.ylabel("Ratio")
  plt.savefig("!{fplot}", dpi=300)
  '''
}
