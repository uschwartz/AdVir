// load modules
// prepare for DANPOS
include{pool} from '../modules/prepareDANPOS'
// DANPOS run
include{danpos_mono; danpos_sub} from '../modules/DANPOS'
// convert DANPOS output
include{convert2saf_mono; convert2saf_sub} from '../modules/convert2saf'
// peak centering
include{centerPeaks} from '../modules/centerPeaks'
// get read count per nucleosome
include{featureCounts_mono; featureCounts_sub} from '../modules/featureCounts'
// get nucMACC_scores
include{nucMACC_scores;sub_nucMACC_scores} from '../modules/nucMACC_scores'

// bamEntry module
include{mergeBam_mono; mergeBam_sub} from '../modules/bamEntry'

workflow sub_bamEntry{

    take:
    bamEntry_mono

    main:
    mergeBam_mono(bamEntry_mono.groupTuple())

    emit:
    mergeBam_mono.out

  }


workflow common_nucMACC{

  take:
  sieve_mono
  samples_conc

  main:
  //get lowest MNase digest
  samples_conc.map{conc,sample -> conc}.min().set{min_conc}

  // get sample with lowest MNase digest
  min_conc.join(samples_conc)
  .map{conc,sample -> sample}
  .set{min_conc_sample}

  // monoNucs
  pool(sieve_mono.map{name,bam -> file(bam)}.collect())
  danpos_mono(sieve_mono.mix(pool.out[0]), pool.out[1])
  convert2saf_mono(danpos_mono.out[1].join(pool.out[0]))
  centerPeaks(convert2saf_mono.out[0], convert2saf_mono.out[1])
  featureCounts_mono(centerPeaks.out[1], sieve_mono.map{name,bam -> file(bam)}.collect())

  //nucMACC scores
  nucMACC_scores(featureCounts_mono.out[0], Channel.fromPath(params.csvInput),featureCounts_mono.out[1])

}
