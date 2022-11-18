process sieve{
  label 'big'
  publishDir "${params.outDir}/QC/04_ALIGNMENT_FILTERING/${anti}", mode: 'copy', pattern: "*_FiltLog.txt"


  input:
  tuple val(sampleID), file(bam), file(bai), val(anti)

  output:
  tuple val(sampleID), file("*_FiltLog.txt")
  tuple val(sampleID), file("*_filt.bam"), val(anti)

  script:
  """
  alignmentSieve -b $bam \
  -o ${sampleID}"_filt.bam" \
  -p $task.cpus \
  --filterMetrics  ${sampleID}"_FiltLog.txt" \
  --minFragmentLength $params.minLen \
  --maxFragmentLength $params.maxLen \
  """
}
