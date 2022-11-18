AnalysisDir=~/R008_ChIP_AdVIr/NFrun_220705_both

cd $AnalysisDir


multiBigwigSummary bins -b RUN/02_VisualizeBAM/Flag/*.bw \
 RUN/02_VisualizeBAM/H3K27ac/*.bw RUN/02_VisualizeBAM/Input/*.bw \
-o downstream/02_Exploratory/results.npz \
--smartLabels --binSize 50 -r AdVir_del:0:34062 -p 12 \
--outRawCounts downstream/02_Exploratory/counts_multiBW.txt


plotPCA --transpose --corData downstream/02_Exploratory/results.npz \
-o downstream/02_Exploratory/PCA.pdf \
 --outFileNameData downstream/02_Exploratory/PCA.tab

plotCorrelation --corData downstream/02_Exploratory/results.npz \
 --corMethod pearson --whatToPlot heatmap \
  -o downstream/02_Exploratory/CorHeatmap.pdf \
    --outFileCorMatrix downstream/02_Exploratory/CorHeatmap.tab
