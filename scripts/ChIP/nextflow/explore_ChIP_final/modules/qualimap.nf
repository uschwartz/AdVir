process qualimap{
  label 'mid'
  publishDir "${params.outDir}/QC/03_qualimap/${anti}", mode: 'copy'


  input:
  tuple val(sampleID), file(bam), file(bam_idx), val(anti)

  output:
  file "${sampleID}"

  script:
  """
  export JAVA_HOME=`/usr/libexec/java_home -v 1.8`
  qualimap bamqc --java-mem-size=16G -bam $bam -c -outdir ${sampleID}
  """

}
