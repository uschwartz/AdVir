process idxstats{
  publishDir "${params.outDir}/QC/02_IDXSTATS/${anti}", mode: 'copy'


  input:
  tuple val(sampleID), file(bam), file(bai), val(anti)

  output:
  file "*_idxStats.txt"


  script:
  """
  samtools idxstats $bam >${sampleID}'_idxStats.txt'
  """
}
