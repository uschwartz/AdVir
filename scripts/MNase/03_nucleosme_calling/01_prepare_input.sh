# set directory where aligned MNase-seq data is stored
alignMNaseDir="~/MNase/alignment"



cd $alignMNaseDir/bed
mkdir -p nucleR

for bed in *hi*; do

# take start and ends of fragment and convert it to 6-col bed format
name=${bed:0:${#bed}-10}
awk '{print $1, $2, $6, $7, $8, $9}' $name"_bedpe.bed" | tr ' ' '\t' >nucleR/$name".bed"


done

