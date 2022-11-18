AnalysisDir=~/R008_ChIP_AdVIr

cd $AnalysisDir


nextflow run nextflow/explore_ChIP_final \
--csvInput 'additionalFiles/fastqIN_220705_both.csv' \
--outDir 'NFrun_220705_both' \
--genomeIdx $AnalysisDir'/data/Annotation/bwt2idx/hg19_AdVir_del' \
-w 'work_220705_both' -resume
