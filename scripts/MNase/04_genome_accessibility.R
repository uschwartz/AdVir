library(GenomicRanges)
library(rtracklayer)

## set working directory
inputDataDir="~/MNase/pVII/quantile_norm_low/bigWig"
setwd(inputDataDir)


## load pVII low MNase profiles

low.0.r1<-import.bw("BorV2_T0_low.bw")

low.05.r1<-import.bw("BorV1_T0.5_low.bw")
low.05.r2<-import.bw("BorV2_T0.5_low.bw")

######## exchange r1 and r2 ######################
low.1.r1<-import.bw("BorV1_T1_low.bw")
low.1.r2<-import.bw("BorV2_T1_low.bw")

low.4.r1<-import.bw("BorV1_T4_low.bw")
low.4.r2<-import.bw("BorV2_T4_low.bw")


######################### get average signal across replicates #######
# convert Granges to Rle
gr2Rle<-function(gr){
    rle.obj<-coverage(gr, weight = values(gr)$score)$AdVir_del
    return(rle.obj)
}

#function Rle to  Granges
rle2gr<-function(rle){
    gr<-GRanges(seqnames = "AdVir_del",ranges=IRanges(start(rle),
    end(rle)), score=runValue(rle))
    return(gr)
}

R1.0.rle<-gr2Rle(low.0.r1)

R1.0.5.rle<-gr2Rle(low.05.r1)
R2.0.5.rle<-gr2Rle(low.05.r2)

R1.1.rle<-gr2Rle(low.1.r1)
R2.1.rle<-gr2Rle(low.1.r2)

R1.4.rle<-gr2Rle(low.4.r1)
R2.4.rle<-gr2Rle(low.4.r2)

####### get average
av.0.5.rle<-(R1.0.5.rle+R2.0.5.rle)/2
av.1.rle<-(R1.1.rle+R2.1.rle)/2
av.4.rle<-(R1.4.rle+R2.4.rle)/2


############################# Calculate accessibility score ####################################
##### define window size and sliding
win.size<-300
win.slide<-20

total<-34063
starts<-seq(1,total-win.size, win.slide)
ends<-seq(win.size,total, win.slide)

#calculate average signal in window
win.0<-aggregate(R1.0.rle,IRanges(starts,end = ends),FUN=mean)
win.0.5<-aggregate(R1.0.5.rle,IRanges(starts,end = ends),FUN=mean)
win.1<-aggregate(R1.1.rle,IRanges(starts,end = ends),FUN=mean)
win.4<-aggregate(R1.4.rle,IRanges(starts,end = ends),FUN=mean)

#merge to table
df.win<-(cbind(win.0, win.0.5, win.1, win.4))

## get slope (=accessibility score)
getSlope<-function(x){
  time.points<-c(0,0.5,1,4)
  fit<-lm(x~time.points)
  slope<-fit$coefficients[2]
}

# apply to all windows
slopes<-apply(df.win, 1, getSlope)

# calculate goodness of fit
getRsqr<-function(x){
  time.points<-c(0,0.5,1,4)
  fit<-lm(x~time.points)
  R2<-summary(fit)$adj.r.squared
}
Rsqr<-apply(df.win, 1, getRsqr)
regr.access<-data.frame(df.win, slopes, Rsqr)
save(regr.access, file="regr.access.rda")

#getz the center of the window
shift<-(win.size-win.slide)/2
sub.df<-data.frame(start=starts+shift, end=ends-shift, regr.access$slopes)

#convert accessibility score to profile (bedgraph)
bedGraph<-data.frame(seqnames=rep("AdVir_del",nrow(sub.df) ) ,sub.df )
write.table(bedGraph,"regr_analysis.bedGraph" ,col.names = F, row.names = F, 
            sep="\t", quote=F)






