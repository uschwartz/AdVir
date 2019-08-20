#### Configure paths

## set directory where annotation is stored
annotationDir="~/annotation"

# set directory where raw MNase-seq data is stored
rawMNaseDir="~/data/MNase"

# set directory of alignment output
outMNaseDir="~/MNase/alignment"


########## align sequencing reads

cd $rawMNaseDir

# run for each fastq pair
for entry in *1_sequence.txt.gz; do

echo 'mapping for sample' $entry

# get main ID
sampleID=$(echo $entry| cut -d'_' -f 2-5)
sampleName=$(echo $entry| cut -d'_' -f 1-8)


# map to genome
bowtie2 -x $annotationDir/Bw_idx/AdVir_del_  -1 $sampleName"_1_sequence.txt.gz"  -2 $sampleName"_2_sequence.txt.gz"   -S $outMNaseDir/$sampleID".sam"  -p 12 --local -t --very-sensitive-local --no-discordant

# filter MAPQC < 30 pairs and convert to bam
samtools view -bS -f 2 -q 30 $outMNaseDir/$sampleID".sam" >$outMNaseDir/$sampleID"_q30.bam";


# get only fragments mapping to AdVir_del
samtools sort  -T $outMNaseDir/$sampleID"_sorted" -o $outMNaseDir/$sampleID"_sorted.bam" $outMNaseDir/$sampleID"_q30.bam"
samtools index $outMNaseDir/$sampleID"_sorted.bam";

samtools view -b $outMNaseDir/$sampleID"_sorted.bam" AdVir_del >$outMNaseDir/$sampleID"_AdVir.bam"


# remove intermediate files
rm $outMNaseDir/$sampleID".sam"
rm $outMNaseDir/$sampleID"_q30.bam"

done


