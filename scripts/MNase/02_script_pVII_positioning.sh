#### Configure paths

## set directory where annotation is stored
annotationDir="~/annotation"

# set directory where aligned MNase-seq data is stored
alignMNaseDir="~/MNase/alignment"

# set directory where aligned MNase-seq data is stored
outDANPOSdir="~/MNase/pVII"


########## align sequencing reads

cd alignMNaseDir
mkdir -p bed
mkdir -p bed/DANPOS_input
mkdir -p bed/DANPOS_input
mkdir -p $outDANPOSdir


# run for each fastq pair
for entry in *"_AdVir.bam"; do

echo 'processing for sample' $entry

# get main ID
sampleID=$(echo $entry| cut -d'_' -f 1-3)

# sort by name
samtools sort -n $entry $sampleID"_sort"
rm $entry

# make size selection; all fragments shorter than 140bp
bamToBed -bedpe -i $sampleID"_sort.bam" >bed/$sampleID"_bedpe.bed"
awk '$6-$2 < 140' bed/$sampleID"_bedpe.bed" >bed/$sampleID"_bedpe_size.bed"

# convert back to bed
bedtools bedpetobam -i bed/$sampleID"_bedpe_size.bed"  -g $annotationDir/ChromAdVir_del_info.txt >bed/$sampleID"_conversion.bam"


bedtools bamtobed -i bed/$sampleID"_conversion.bam" >bed/DANPOS_input/$sampleID"_re-conversion.bed"

# remove intermediate files
rm bed/$sampleID"_bedpe_size.bed"
rm bed/$sampleID"_conversion.bam"

## get pVII positions and profile at AdVir_del

python danpos.py dpos bed/DANPOS_input/$sampleID"_re-conversion.bed" -m 1 -c 20000 -u 1e-5 -a 1 -jd 50 -z 5 -p 0.001 --extend 35 --mifrsz 0 -o $outDANPOSdir/$sampleID >$outDANPOSdir/$sampleID/running_info.txt

done


#### quantile normalization

##################
## high MNase
##################
#create directory
mkdir -p $outDANPOSdir/quantile_norm_high
# set reference profile to high MNase 0 hpi Replica 1
reference="BorV1_T0_hi"

python danpos.py wig2wiq $annotationDir/ChromAdVir_del_info.txt $outDANPOSdir/$reference/pooled/*.Fnor.smooth.wig --out_dir $outDANPOSdir/quantile_norm_high

## normalize to reference
while read name; do
echo $name

python danpos.py wiq $annotationDir/ChromAdVir_del_info.txt $outDANPOSdir/$name/pooled/*.Fnor.smooth.wig --reference  $outDANPOSdir/quantile_norm_high/$reference*".wiq"  --step 1 --out_dir $outDANPOSdir/quantile_norm_high/ --rformat wiq --rsorted 1

done < $annotationDir/quantile_norm/cond_high.sh


##################
## low MNase
##################
mkdir -p $outDANPOSdir/quantile_norm_low

# set reference profile to low MNase 0 hpi Replica 1
reference="BorV2_T0_low"

python danpos.py wig2wiq $annotationDir/ChromAdVir_del_info.txt $outDANPOSdir/$reference/pooled/*.Fnor.smooth.wig --out_dir $outDANPOSdir/quantile_norm_low

## normalize to reference
while read name; do
echo $name

python danpos.py wiq $annotationDir/ChromAdVir_del_info.txt $outDANPOSdir/$name/pooled/*.Fnor.smooth.wig --reference  $outDANPOSdir/quantile_norm_low/$reference*".wiq"  --step 1 --out_dir $outDANPOSdir/quantile_norm_low/ --rformat wiq --rsorted 1

done < $annotationDir/quantile_norm/cond_low.sh



########### get pVII positions from normalized profiles ######

##################
## high MNase
##################
cd $outDANPOSdir/quantile_norm_high/

mkdir -p bigWig

for wig in *qnor.wig;do

name=$(echo $wig| cut -d'_' -f 1-3)

mkdir $name

python danpos.py dpos $wig -jd 50 -z 5 -q 10 -a 1 -o $name >$name/running_info.txt

wigToBigWig $name/pooled/*wig $annotationDir/ChromAdVir_del_info.txt "bigWig/"$name".bw"

done

##################
## low MNase
##################
cd $outDANPOSdir/quantile_norm_low/

mkdir -p bigWig

for wig in *qnor.wig;do

name=$(echo $wig| cut -d'_' -f 1-3)

mkdir $name

python danpos.py dpos $wig -jd 50 -z 5 -q 10 -a 1 -o $name >$name/running_info.txt

wigToBigWig $name/pooled/*wig $annotationDir/ChromAdVir_del_info.txt "bigWig/"$name".bw"

done
