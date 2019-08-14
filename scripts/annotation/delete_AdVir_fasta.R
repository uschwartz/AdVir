library(seqinr)

# Human adenovirus C serotype 5, complete genome was downloaded from NCBI; Accession AY339865.1
AdVir<-read.fasta("AdVir_C_wt_genome.fasta", forceDNAtolower = F)

# The virus strain used has a deletion in E3 region: 28593-30464bp
# write deleted version to output
write.fasta(AdVir[[1]][-(28593:30464)], names="AdVir_del", file.out="AdVir_del.fa")

