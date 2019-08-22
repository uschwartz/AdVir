setwd("~/RNA/alignment/")
library(Gviz)
library(GenomicRanges)
library(rtracklayer)
options(ucscChromosomeNames=FALSE)
gtrack <- GenomeAxisTrack()

## get gene annotation
library(GenomicFeatures)
AdVir.TxDb<-makeTxDbFromGFF("~/annotation/AdVir_del_final.gtf",format="gtf")
grTrackGene<-GeneRegionTrack(AdVir.TxDb, transcriptAnnotation="gene", 
                             collapseTranscripts = TRUE)

################ get Data ###
#take only AdVir choromosome
AdVir<-GRanges(seqnames = "AdVir_del", ranges = IRanges(1,34062))

# get profiles of uniquely mapped reads
bw.files<-list.files("bigwigsUnique/")
names<-sapply(strsplit(bw.files,split = "_Signal"), function(x) x[1])

list.ExpRLe<-list()
for(i in bw.files){
  name<-unlist(strsplit(i,split = "_Signal"))[1]
  list.ExpRLe[[name]]<-import.bw(con=paste0("bigwigsUnique/",i),
                                 selection=BigWigSelection(AdVir),as="RleList")$AdVir_del
}

for(i in 2:length(list.ExpRLe)){
  if(i==2){
    summedRLE<-list.ExpRLe[[1]]+list.ExpRLe[[i]]
  } else {
    summedRLE<-summedRLE+list.ExpRLe[[i]]
  }
}

######## convert to GRanges obj
GRangesExp<-GRanges(seqnames = "AdVir_del",ranges=IRanges(start(summedRLE),
                                                          end(summedRLE)))
values(GRangesExp)<-sapply(list.ExpRLe,function(x) aggregate(x ,
                                              ranges(GRangesExp),FUN=mean))

######################load data into Gviz
#colors
brew<-RColorBrewer::brewer.pal(5, "Set3")

dTrack <- DataTrack(GRangesExp, name = "normalized expression", 
                    col=brew)

# log-transform plot
logTrans<-function(x){
  x<-x+1
  x[x < 1] <- 1 
  log10(x)
}


# plot data

plotTracks(list(dTrack,grTrackGene,gtrack), sizes = c(6,1, 1.5),
           groups = rep(c("0.5hpi", "0hpi","1hpi", "2hpi", "4hpi"),each=2),
           type = c("a"), transformation = logTrans)




#################### Replicates #########################

#create GRanges obj split by replicate
GRangesExp.r1<-GRangesExp
values(GRangesExp.r1)<-sapply(list.ExpRLe[grep("_r1", names(list.ExpRLe))],
                              function(x) aggregate(x ,ranges(GRangesExp.r1),FUN=mean))

GRangesExp.r2<-GRangesExp
values(GRangesExp.r2)<-sapply(list.ExpRLe[grep("_r2", names(list.ExpRLe))],
                              function(x) aggregate(x ,ranges(GRangesExp.r2),FUN=mean))

#load data into Gviz
dTrack.r1 <- DataTrack(GRangesExp.r1, name = "Replica 1", 
                       col=brew,
                       groups = rep(c("0.5hpi", "0hpi","1hpi", "2hpi", "4hpi")))
dTrack.r2 <- DataTrack(GRangesExp.r2, name = "Replica 2", 
                       col=brew,
                       groups = rep(c("0.5hpi", "0hpi","1hpi", "2hpi", "4hpi")))


#plot data

plotTracks(list(dTrack.r1,dTrack.r2,grTrackGene,gtrack), 
           sizes = c(4,4,1,1.5),
           type = c("a"), transformation = logTrans)


