process centerPeaks{
  container 'uschwartz/r_nucmacc:v3.1'
  publishDir "${params.outDir}/RUN/03_NUCS_POSITIONS", mode: 'copy', pattern: "*_nucPositions_centered.bed"

  input:
  file(bed)
  file(saf)

  output:
  file("*_nucPositions_centered.bed")
  file("*_nucPositions_centered.saf")

  script:
  """
  centerNucs.R
  """
}
