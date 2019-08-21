library(nucleR)
library(rtracklayer)

## set colors
library(RColorBrewer)
purp<-brewer.pal(5,"Purples")

## set working directory
inputDataDir="~/MNase/alignment/bed/nucleR"
setwd(inputDataDir)
dir.create("plots")

## load data from bed
## high MNase
r1.1h.bed<-import.bed(con="BorV1_T1_hi_mono.bed")
r1.2h.bed<-import.bed(con="BorV1_T2_hi_mono.bed")
r1.4h.bed<-import.bed(con="BorV1_T4_hi_mono.bed")

ir.1h.r1<- as(r1.1h.bed, "RangedData")
ir.2h.r1<- as(r1.2h.bed, "RangedData")
ir.4h.r1<- as(r1.4h.bed, "RangedData")

# plot fragment length
pdf("plots/time_course_fragment_length.pdf")
    plot(density(width(ir.1h.r1)), col=purp[2], lwd=3, main="nucleosome assembly",
         xlab="fragment length", ylim=c(0,0.017))
    lines(density(width(ir.2h.r1)), col=purp[4], lwd=3)
    lines(density(width(ir.4h.r1)), col=purp[5], lwd=3)
    abline(v=c(137,157))
    legend("topright", legend = c("1 hpi", "2 hpi", "4 hpi"), bty="n",
           col=pup[c(2,4,5)], lwd=3)
dev.off()

#extract fragments with length of 137-157 bp
nucs.4h.r1<-ir.4h.r1[(width(ir.4h.r1)>=137 & width(ir.4h.r1)<=157), ]
nucs.2h.r1<-ir.2h.r1[(width(ir.2h.r1)>=137 & width(ir.2h.r1)<=157), ]
nucs.1h.r1<-ir.1h.r1[(width(ir.1h.r1)>=137 & width(ir.1h.r1)<=157), ]


#calculate proportion of fragments with 139-159bp length
quant.nuc.4h.r1<-length(ranges(nucs.4h.r1)$AdVir)
quant.nuc.2h.r1<-length(ranges(nucs.2h.r1)$AdVir)
quant.nuc.1h.r1<-length(ranges(nucs.1h.r1)$AdVir)

quant.all.4h.r1<-length(ranges(ir.4h.r1)$AdVir)
quant.all.2h.r1<-length(ranges(ir.2h.r1)$AdVir)
quant.all.1h.r1<-length(ranges(ir.1h.r1)$AdVir)

prop.nucs.4h.r1<-quant.nuc.4h.r1/quant.all.4h.r1
prop.nucs.2h.r1<-quant.nuc.2h.r1/quant.all.2h.r1
prop.nucs.1h.r1<-quant.nuc.1h.r1/quant.all.1h.r1

#process Reads
reads_pair.r1.4h<- processReads(nucs.4h.r1, type="paired")
reads_pair.r1.2h<- processReads(nucs.2h.r1, type="paired")
reads_pair.r1.1h<- processReads(nucs.1h.r1, type="paired")



#calculate coverage and normalize to rpm
cover.4h.r1 <- coverage.rpm(reads_pair.r1.4h, scale=1000000*prop.nucs.4h.r1)
cover.2h.r1 <- coverage.rpm(reads_pair.r1.2h, scale=1000000*prop.nucs.2h.r1)
cover.1h.r1 <- coverage.rpm(reads_pair.r1.1h, scale=1000000*prop.nucs.1h.r1)



#Overlay
palette(brewer.pal(8,"Set1"))


nuc.4h = as.vector(cover.4h.r1[[1]])
nuc.2h = as.vector(cover.2h.r1[[1]])
nuc.1h =as.vector(cover.1h.r1[[1]])

htseq_raw.4h.r1<-c(nuc.4h,rep(0,34062-length(nuc.4h)))
htseq_raw.2h.r1<-c(nuc.2h,rep(0,34062-length(nuc.2h)))
htseq_raw.1h.r1<-c(nuc.1h,rep(0,34062-length(nuc.1h)))

htseq_fft.4h.r1 = filterFFT(htseq_raw.4h.r1,pcKeepComp=0.007, showPowerSpec = T)
htseq_fft.2h.r1 = filterFFT(htseq_raw.2h.r1, pcKeepComp=0.007,showPowerSpec = T)
htseq_fft.1h.r1 = filterFFT(htseq_raw.1h.r1, pcKeepComp=0.007,showPowerSpec = T)


## Visualization

palette(brewer.pal(8, "Set1"))

pdf(file="plots/raw_data_Difference_Nois_reduction.pdf",width=11.6, 
        height=8.2 )
    par(mfrow=c(3,1))
    plot(htseq_raw.1h.r1, type="l", col=1, lwd=3, ylab="coverage",
         xlab="position", main="Replica2 1hpi",
         ylim=range(htseq_raw.4h.r1,htseq_raw.2h.r1), bty="n" )
    lines(htseq_fft.1h.r1, type="l", col="black")
    legend("topleft", legend=c("noise filtered", "raw data"), col=c( "black", 1),
           lty=1 , bty="n")
    plot(htseq_raw.2h.r1, type="l", col=2, lwd=3, ylab="coverage",
         xlab="position", main="Replica2 2hpi",
         ylim=range(htseq_raw.4h.r1,htseq_raw.2h.r1), bty="n" )
    lines(htseq_fft.2h.r1, type="l", col="black")
    legend("topleft", legend=c("noise filtered", "raw data"), col=c( "black", 2),
           lty=1 , bty="n")
    plot(htseq_raw.4h.r1, type="l", col=3, lwd=3, ylab="coverage",
         xlab="position", main="Replica2 4hpi",
         ylim=range(htseq_raw.4h.r1,htseq_raw.2h.r1), bty="n" )
    lines(htseq_fft.4h.r1, type="l", col="black")
    legend("topleft", legend=c("noise filtered", "raw data"), col=c( "black", 3),
           lty=1 , bty="n")
dev.off()




## peak calling
peaks.2h.r1=peakDetection(htseq_fft.2h.r1, score=TRUE, width=147)
peaks.4h.r1=peakDetection(htseq_fft.4h.r1, score=TRUE, width=147)


pdf(file="plots/peak_calling_ZOOM.pdf",width=11.6,  height=8.2 )
    par( mfrow=c(3,1))
    plotPeaks(peaks.4h.r1, htseq_fft.4h.r1, threshold=400, ylim=c(0, 2500),
              xlim=c(15000,25000),bty="n", main="Replica 4hpi - Peak calling")
    plotPeaks(peaks.4h.r1, htseq_fft.4h.r1, threshold=400, ylim=c(0, 2500),
              xlim=c(18000,22000),bty="n", main="Replica 4hpi - Peak calling")
    plotPeaks(peaks.4h.r1, htseq_fft.4h.r1, threshold=400, ylim=c(0, 1500),
              xlim=c(18000,19000),bty="n", main="Replica 4hpi - Peak calling")
dev.off()


peaks.4h.r1


### 2hpi
## calculate signal under peaks and compare to background at 1hpi
val.2h<-aggregate(as(htseq_fft.2h.r1,"Rle"),FUN=mean, start=start(peaks.2h.r1),end=end(peaks.2h.r1))
input.2h<-aggregate(as(htseq_fft.1h.r1,"Rle"),FUN=mean, start=start(peaks.2h.r1),end=end(peaks.2h.r1))

# peaks with at least 5 fold change increase in signal
valid.2h.5<-val.2h/input.2h>5
peaks.2h.r1.filt5<-peaks.2h.r1[valid.2h.5,]


### 4hpi
## calculate signal under peaks and compare to background at 1hpi
val.4h<-aggregate(as(htseq_fft.4h.r1,"Rle"),FUN=mean, start=start(peaks.4h.r1),end=end(peaks.4h.r1))
input.4h<-aggregate(as(htseq_fft.1h.r1,"Rle"),FUN=mean, start=start(peaks.4h.r1),end=end(peaks.4h.r1))

# peaks with at least 5 fold change increase in signal
valid.4h.5<-val.4h/input.4h>5
peaks.4h.r1.filt5<-peaks.4h.r1[valid.4h.5,]


#write bed files
names(peaks.4h.r1)<-"AdVir_del"
names(peaks.2h.r1)<-"AdVir_del"

#all peaks detected
rtracklayer::export.bed(peaks.2h.r1, con="peaks/NucPeaks_replica2_2h.bed")
rtracklayer::export.bed(peaks.4h.r1, con="peaks/NucPeaks_replica2_4h.bed")

#filtered peaks with at least 5 fold change
rtracklayer::export.bed(peaks.2h.r1.filt5, con="peaks/NucPeaks_replica2_2h_fc5.bed")
rtracklayer::export.bed(peaks.4h.r1.filt5, con="peaks/NucPeaks_replica2_4h_fc5.bed")


#export profile of nucleosome sized fragments
gr.1h.r1<-GRanges(seqnames="AdVir_del", ranges=IRanges(1:34062,width=1), score=htseq_fft.1h.r1)
gr.2h.r1<-GRanges(seqnames="AdVir_del", ranges=IRanges(1:34062,width=1), score=htseq_fft.2h.r1)
gr.4h.r1<-GRanges(seqnames="AdVir_del", ranges=IRanges(1:34062,width=1), score=htseq_fft.4h.r1)
rtracklayer::export.wig(gr.1h.r1, con="Replica2_1_hpi_nucs.wig")
rtracklayer::export.wig(gr.2h.r1, con="Replica2_2_hpi_nucs.wig")
rtracklayer::export.wig(gr.4h.r1, con="Replica2_4_hpi_nucs.wig")



