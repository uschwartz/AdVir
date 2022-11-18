process multiqc{
  publishDir "${params.outDir}/QC/multiqc/", mode: 'copy'

  input:
  file('*')
  output:
  file "*multiqc_report.html"
  file "*_data"

  script:
  """
  multiqc .
  """
}
