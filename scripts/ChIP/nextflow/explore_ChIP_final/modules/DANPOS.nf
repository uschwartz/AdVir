process danpos{
  container 'uschwartz/danpos'
  memory { params.genomeSize > 200000000 ? 50.GB : 16.GB }
  publishDir "${params.outDir}/RUN/03_DANPOS_PROFILE/${anti}", mode: 'copy',
  pattern: "*_profile.bw"


  input:
  tuple val(sampleID), file(bam), val(anti)
  file(chrSizes)

  output:
  file("*_profile.bw")
  tuple val(sampleID), file("result/pooled/*.xls")

  script:
  """
  danpos.py dpos $bam -m 1 --extend 70 -c $params.genomeSize \
  -u 0 -z 20 -e 1 \
  --distance 75 --width 10  > $sampleID"_DANPOS_stats.txt"
  wigToBigWig result/pooled/*.wig -clip $chrSizes $sampleID"_profile.bw"
  """
}
