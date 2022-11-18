process bam2bw{
  label 'mid'
  publishDir "${params.outDir}/RUN/02_VisualizeBAM/${anti}", mode: 'copy'


  input:
  tuple val(sampleID), file(bam), file(bam_idx), val(anti)

  output:
  tuple val(sampleID), file("*.bw")

  script:
  """
  bamCoverage -b $bam \
   -o ${sampleID}".bw" \
   -p $task.cpus \
   --normalizeUsing 'CPM' \
   --extendReads
  """

}
