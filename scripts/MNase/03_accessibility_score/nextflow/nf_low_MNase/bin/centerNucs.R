#!/usr/bin/env Rscript

saf.file<-list.files(pattern = "*.saf")


saf<-read.delim(saf.file)

saf$Start<-saf$Start+39
saf$End<-saf$End-39

write.table(saf, "pooled_nucPositions_centered.saf",
            row.names = F, col.names = T, quote=F,
            sep="\t")


bed.file<-list.files(pattern = "*.bed")

bed<-read.delim(bed.file, header = F)

bed$V2<-bed$V2+39
bed$V3<-bed$V3-39

write.table(bed, "pooled_nucPositions_centered.bed",
            row.names = F, col.names = F, quote=F,
            sep="\t")
