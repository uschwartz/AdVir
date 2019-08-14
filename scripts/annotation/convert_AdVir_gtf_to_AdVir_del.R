gtf<-read.delim("AdVir_del_flattend.gtf",header=F)

#change chr
gtf$V1<- "AdVir_del"

#delet parts
start_del<-28593
end_del<-30464
diff<-end_del-start_del

## start and end in deletion -> shift back by length deletion
del.idx<-which(gtf$V4>start_del & gtf$V5<end_del )

gtf<-gtf[-del.idx,]



###########
### shift exon ends in deletion
end.idx<-which(gtf$V4<start_del & gtf$V5<end_del & gtf$V5>start_del)
gtf$V5[end.idx]<-start_del

### shift exon starts in deletion
start.idx<-which(gtf$V4>start_del & gtf$V4<end_del & gtf$V5>end_del)
gtf$V4[start.idx]<-start_del
gtf$V5[start.idx]<-gtf$V5[start.idx]-diff


## start and end after deletion -> shift back by length deletion
shift.idx<-which(gtf$V4>end_del)

gtf$V4[shift.idx]<-gtf$V4[shift.idx] - diff
gtf$V5[shift.idx]<-gtf$V5[shift.idx] - diff


##############
write.table(gtf, file="AdVir_del_final.gtf", quote=F, col.names = F,
            row.names = F, sep="\t")
