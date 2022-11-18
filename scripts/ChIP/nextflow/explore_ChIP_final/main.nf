#!/usr/bin/env nextflow

 /*
 ===============================================================================
             nextflow based ChIP-seq explore analysis pipeline
 ===============================================================================
Authors:
Uwe Schwartz <uwe.schwartz@ur.de>
 -------------------------------------------------------------------------------
 */

nextflow.enable.dsl = 2

 //                           show settings
 if (!params.help) {
         include{settings} from './modules/setting'
         settings()
 }



 //                      workflow
// read csv file
// forward reads
Channel.fromPath(params.csvInput)
        .splitCsv(header:true)
        .map{ row -> tuple(row.name,file(row.path_fwdReads), row.condition, row.antibody)}
        .set{samples_fwd_ch}

// reverse reads
Channel.fromPath(params.csvInput)
      .splitCsv(header:true)
      .map{ row -> tuple(row.name,file(row.path_revReads),row.condition, row.antibody)}
      .set{samples_rev_ch}



//Channel for fastqc
samples_fwd_ch.mix(samples_rev_ch).set{sampleSingle_ch}
samples_fwd_ch.map{
        name,path_fwdReads,condition,antibody -> tuple(name,file(path_fwdReads))
        }.join(samples_rev_ch).set{samplePair_ch}

////////////////////////////////////////////////////////////

// load modules
//fastqc
include{fastqc} from './modules/raw_qc'
//multiqc
include{multiqc} from './modules/multiqc'
//alignment to ref genome
include{alignment} from './modules/align'
//idxStats
include{idxstats} from './modules/idxstats'
//bam2bw
include{bam2bw} from './modules/bam2bw'
//qualimap
include{qualimap} from './modules/qualimap'
// filtering sizes using alignmentSieve
include{sieve} from './modules/alignmentsieve'
// prepare for DANPOS
include{chrsize} from './modules/prepareDANPOS'
// DANPOS run
include{danpos} from './modules/DANPOS'


workflow{
  fastqc(sampleSingle_ch)
  alignment(samplePair_ch)
  idxstats(alignment.out[1])

  bam2bw(alignment.out[1])

  qualimap(alignment.out[1])
  //Fragment sizes
  // get mono-nucleosomal fragments
  sieve(alignment.out[1])
  //get chrom_Sizes
  chrsize(sieve.out[1].map{name,bam,abti -> file(bam)}.collectFile(sort:true)
  .first())
  //get nucleosome profiles
  danpos(sieve.out[1], chrsize.out)


  multiqc(fastqc.out[0].mix(alignment.out[0]).mix(qualimap.out).collect())

}
