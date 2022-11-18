AnalysisDir=~/R008_ChIP_AdVIr/data/lowMNase/bam_sorted

cd $AnalysisDir

for entry in *sort.bam; do

    mainID=$(echo $entry | cut -d'_' -f -7 )
    echo $mainID
    
    samtools sort $entry >$mainID".bam"
    samtools index $mainID".bam"
    
    alignmentSieve -b $mainID".bam" \
     -p 5 --maxFragmentLength 140 \
      -o .."/"$mainID"_sub.bam"
      
done
