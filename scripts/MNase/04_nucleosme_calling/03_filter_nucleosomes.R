setwd("~/R008_ChIP_AdVIr/NFrun_220705_both/downstream/05_Nucs/")
library(nucleR)
library(rtracklayer)

## set colors
library(RColorBrewer)
purp<-brewer.pal(5,"Purples")

### read flag H3.3
chip.table<-read.delim("../02_Exploratory/counts_multiBW.txt")

flag.t4<-grep("U2OS33_T4h_Flag",colnames(chip.table), value = T)
flag.t0<-grep("U2OS33_T0h_Flag",colnames(chip.table), value = T)

dist.t4<-apply(chip.table[,flag.t4],1,mean)

enrch.regions<-chip.table[dist.t4>2,1:3]

enrich.gr<-GRanges(seqnames =enrch.regions$X..chr.,
        ranges=IRanges(start=enrch.regions$X.start.,
                       end =enrch.regions$X.end. ))
#merge consecutive elements
enrch.flatten<-reduce(enrich.gr)


### get nucleosomes
nucs.rep1<-import.bed("~/MNase/alignment/bed/nucleR/NucPeaks_replica1_4h.bed")
nucs.rep2<-import.bed("~/MNase/alignment/bed/nucleR/NucPeaks_replica2_4h.bed")


ovrl.r2<-findOverlaps(nucs.rep2,enrch.flatten, type = "within")
export.bed(nucs.rep2[ovrl.r2@from], "nucs_rep2_flt.bed")

ovrl.r1<-findOverlaps(nucs.rep1,enrch.flatten, type = "within")
export.bed(nucs.rep1[ovrl.r1@from], "nucs_rep1_flt.bed")

export.bed(enrch.flatten, "H3.3_enriched-region.bed")
