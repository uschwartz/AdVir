AnalysisDir=~/R008_ChIP_AdVIr


##run script CHANGED FEATURECOUNTS

cd $AnalysisDir

nextflow run  $AnalysisDir'/nextflow/nf_low_MNase' \
--analysis 'nucMACC'  --bamEntry \
--csvInput $AnalysisDir'/addData/sample_sheet_lowMNase.csv' \
--outDir $AnalysisDir'/lowMNase_nf' \
--genomeSize 34062 \
--genome $AnalysisDir'/data/Annotation/AdVir_del.fa' \
-w ./work_lowMNase_nf -resume
