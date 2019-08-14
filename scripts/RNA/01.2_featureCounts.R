
args<-commandArgs(TRUE)

library(Rsubread)
bams<-list.files()[grep(".bam",list.files())]
gtf<-args[1]
counts<-featureCounts(files=bams, annot.ext=gtf, nthreads=10,
                      isGTFAnnotationFile=T)

count.table<-(counts[["counts"]])
length.genes<-(counts[["annotation"]]$Length)

#calculate TPMs
count.norm1<-apply(count.table,2,function(x) x/length.genes)
normFac<-1e6/apply(count.norm1,2,sum)

gx<-grep("AdVir_del", counts[["annotation"]]$Chr)
counts[["annotation"]][gx,]

count.tpm.advir<-apply(count.norm1[gx,],1,function(x) x*normFac)


######################### plot data
library(ggplot2)
library(reshape2)
library(ggbeeswarm)

#transform to data frame
hpi<-(rep(c(0.5,0,1,2,4),each=2))
df<-data.frame(count.tpm.advir[,c("E1A", "E1B", "E2B", "E2A", "E3", "E4")],
hpi)

df.metl<-melt(df,id="hpi", variable.name = "gene", value.name = "expr")

#colors
brew<-RColorBrewer::brewer.pal(5, "Set3")

#########################################################
g<-ggplot(df.metl, aes(x=as.factor(hpi),y=log2(expr+1), fill=as.factor(hpi)))
g<-g+geom_bar(stat = "summary", fun.y = "mean", position=position_dodge(), color="black")+
theme_bw()+ylab("log2 ( TPM + 1 )")+xlab("hpi")+facet_wrap(~gene,ncol = 3)+
scale_fill_manual(values=brew[c(2,1,3:5)], name="hpi")+
geom_point(color="black", position = position_dodge(width = .9))

print(g)
