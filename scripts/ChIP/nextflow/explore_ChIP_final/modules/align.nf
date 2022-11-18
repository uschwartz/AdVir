process alignment{
  label 'big'
  publishDir "${params.outDir}/RUN/01_ALIGNMENT/${anti}", mode: 'copy', pattern: '*_aligned.{bam,bam.bai}'
  publishDir "${params.outDir}/QC/02_ALIGNMENT/${anti}", mode: 'copy', pattern: '*_alignment_stats.txt'


  input:
  tuple val(sampleID), file(read1), file(read2), val(cond), val(anti)

  output:
  file "*_alignment_stats.txt"
  tuple val(sampleID), file("*_aligned.bam"), file("*_aligned.bam.bai"), val(anti)

  script:
  """
  bowtie2 -t \
  --threads $task.cpus \
  --very-sensitive-local \
  --no-discordant \
  -x $params.genomeIdx \
  -1 $read1 \
  -2 $read2 \
  2> ${sampleID}_alignment_stats.txt \
  | samtools view -bS -q 30 -f 2 -@ $task.cpus - | samtools sort -@ $task.cpus - > ${sampleID}"_aligned.bam"

  samtools index -b ${sampleID}"_aligned.bam"
  """
}
