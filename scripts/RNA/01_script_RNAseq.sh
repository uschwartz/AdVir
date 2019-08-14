#### Configure paths

## set directory where annotation is stored
annotationDir="~/annotation"

# set directory where raw RNAseq data is stored
rawRNAdir="~/data/RNA"

# set directory of alignment output
outRNAdir="~/RNA/alignment"

# set directory for counting
countRNAdir="~/RNA/count"

# set directory to scripts
scriptsDir="~/scripts/RNA"


############ Create genome index

STAR --runThreadN 12 --runMode genomeGenerate --genomeDir $annotationDir/STAR_idx --genomeFastaFiles $annotationDir/hg19/Sequence/WholeGenomeFasta/genome.fa $annotationDir/AdVir_del.fa --sjdbGTFfile ~$annotationDir/gencode.v25lift37.annotation.gtf --sjdbOverhang  49


########## align sequencing reads

cd $rawRNAdir

# run for each fastq entry
for entry in *.fastq.gz; do

echo 'mapping for sample' $entry

# get main ID
sampleID=$(echo $entry| cut -d'_' -f 2-3)
mkdir $outRNAdir/$sampleID


STAR --outFileNamePrefix $outRNAdir/$sampleID/$sampleID"_" --runThreadN 12 --readFilesCommand gunzip -c --outMultimapperOrder Random --outSAMmultNmax 1 --outSAMtype BAM SortedByCoordinate --outWigType wiggle --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000  --outWigStrand Unstranded --outFilterIntronMotifs RemoveNoncanonical --outSAMstrandField intronMotif --genomeDir $annotationDir/STAR_idx  --readFilesIn $entry


# convert to bigwig

wigToBigWig $outRNAdir/$sampleID/$sampleID"_Signal.UniqueMultiple.str1.out.wig" $annotationDir/ChromAdVir_del_info.txt $outRNAdir/$sampleID/$sampleID"_Signal.UniqueMultiple.out.bw"
wigToBigWig $outRNAdir/$sampleID/$sampleID"_Signal.Unique.str1.out.wig" $annotationDir/ChromAdVir_del_info.txt $outRNAdir/$sampleID/$sampleID"_Signal.Unique.out.bw"

done

########## count sequencing reads
cd $outRNAdir

## copy the bam files to folder
for entry in *"_"*; do
echo $entry

cp $entry/*.bam $countRNAdir/.

done

cd $countRNAdir

#count reads for human and AdVir transcripts

featureCounts -T 12 -t exon -g gene_id -s 0 -a $annotationDir/hg19_genes_and_AdVir_del.gtf -o counts.txt *.bam 2>counts_info.txt


# count reads and normalize to TPM
Rscript $scriptsDir/01.2_featureCounts.R $annotationDir/hg19_genes_and_AdVir_del.gtf
