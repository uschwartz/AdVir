#!/usr/bin/env nextflow

 /*
 ===============================================================================
                      nextflow based nucMACC pipeline
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

 //                       help message
 // Show help message
 if (params.help) {
     include{helpMessage} from './modules/help'
     helpMessage()
     exit 0
 }

 //                      workflow

// read csv file
Channel
    .fromPath(params.csvInput)
    .splitCsv(header:true)
    .map{ row -> tuple(row.Sample_Name,file(row.path_mono))}
    .set{bamEntry_mono}


//read MNase concentration
Channel
      .fromPath(params.csvInput)
      .splitCsv(header:true)
      .map{ row -> tuple(row.MNase_U.toDouble(),row.Sample_Name)}
      .set{samples_conc}


// load workflows
// generate profiles
include{sub_bamEntry; common_nucMACC} from './workflows/nucMACC'


workflow{
  sub_bamEntry(bamEntry_mono)
  common_nucMACC(sub_bamEntry.out[0], samples_conc)

}
