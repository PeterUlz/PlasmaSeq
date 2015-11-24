### R code from vignette source 'exomeCopy.Rnw'

###################################################
### code chunk number 1: exomeCopy.Rnw:51-54
###################################################
library(exomeCopy)
gr <- GRanges(seqname="seq1",IRanges(start=1,end=345))
subdivideGRanges(gr)


###################################################
### code chunk number 2: exomeCopy.Rnw:60-70
###################################################
  plot(0,0,xlim=c(0,500),ylim=c(0,25),type="n",yaxt="n",ylab="", 
       xlab="width of input GRanges object",
       main="Effect of subdivideGRanges")
  abline(v=1:5*100,col="grey")
  for (i in 1:24) {
    gr <- GRanges(seqname="chr1",IRanges(start=1,width=(i*20)))
    sbd.gr <- subdivideGRanges(gr)
    arrows(start(sbd.gr),rep(i,length(sbd.gr)),end(sbd.gr),
           rep(i,length(sbd.gr)),length=.04,angle=90,code=3)
  }


###################################################
### code chunk number 3: exomeCopy.Rnw:76-84
###################################################
target.file <- system.file("extdata","targets.bed",package="exomeCopy")
target.df <- read.delim(target.file,header=FALSE, 
                        col.names=c("seqname","start","end")) 
target <- GRanges(seqname=target.df$seqname, 
                  IRanges(start=target.df$start+1,end=target.df$end))
target
target.sub <- subdivideGRanges(target)
target.sub


###################################################
### code chunk number 4: exomeCopy.Rnw:91-97
###################################################
bam.file <- system.file("extdata","mapping.bam",package="exomeCopy")
scanBamHeader(bam.file)[[1]]$targets
levels(seqnames(target.sub))
rdata <- RangedData(space=seqnames(target.sub),ranges=ranges(target.sub))
rdata[["sample1"]] <- countBamInGRanges(bam.file,target.sub)
rdata


###################################################
### code chunk number 5: exomeCopy.Rnw:104-113
###################################################
reference.file <- system.file("extdata","reference.fa",package="exomeCopy")
target.dnastringset <- scanFa(reference.file,target.sub)
getGCcontent <- function(x) {
  GC.count <- letterFrequency(x,"GC")
  all.count <- letterFrequency(x,"ATGC")
  as.vector(ifelse(all.count==0,NA,GC.count/all.count))
}
rdata[["GC"]] <- getGCcontent(target.dnastringset)
rdata


###################################################
### code chunk number 6: exomeCopy.Rnw:124-127
###################################################
data(exomecounts)
dim(exomecounts)
exomecounts[1:5,1:3]


###################################################
### code chunk number 7: exomeCopy.Rnw:133-136
###################################################
plot(start(exomecounts),exomecounts$HG00551,xlim=c(0.8e6,1.8e6),
     xlab="genomic position",ylab="counts",
     main="HG00551 read counts in exonic ranges")


###################################################
### code chunk number 8: exomeCopy.Rnw:144-149
###################################################
exome.samples <- grep("HG.+",colnames(exomecounts),value=TRUE)
sample.columns <- colnames(exomecounts) %in% exome.samples
C <- as.data.frame(unlist(values(exomecounts)[,sample.columns]))
C.norm <- sweep(C,2,colMeans(C),"/")
exomecounts[["bg"]] <- apply(C.norm,1,median)


###################################################
### code chunk number 9: exomeCopy.Rnw:154-156
###################################################
exomecounts[["GC.sq"]] <- exomecounts$GC^2
exomecounts[["width"]] <- width(exomecounts)


###################################################
### code chunk number 10: exomeCopy.Rnw:201-214
###################################################
set.seed(2)
cnv.type <- c("hom.del","het.del","het.dup","hom.dup")
cnv.probs <- c(.99,.5,.5,.95)
cnv.mult <- c(-1,-1,1,1)
bounds <- IRanges(start=3e6,end=4e6)
for (i in 1:4) {
  samplename <- exome.samples[i]
  contained <- unlist(ranges(exomecounts)) %in% bounds
  O <- exomecounts[[samplename]]
  O[contained] <- O[contained] + (cnv.mult[i] * 
    rbinom(sum(contained),prob=cnv.probs[i],size=O[contained]))
  exomecounts[[paste(samplename,cnv.type[i],sep=".")]] <- O
}


###################################################
### code chunk number 11: exomeCopy.Rnw:218-226
###################################################
par(mfrow=c(2,1),mar=c(5,4,3,2))
plot(start(exomecounts),exomecounts[["HG00731"]],
     xlab="genomic position",ylab="counts",main="Original counts")
abline(v=c(start(bounds),end(bounds)))
plot(start(exomecounts),exomecounts[["HG00731.het.dup"]],
     xlab="genomic position",ylab="counts",
     main="Simulated heterozygous duplication")
abline(v=c(start(bounds),end(bounds)))


###################################################
### code chunk number 12: exomeCopy.Rnw:234-237
###################################################
  fit <- exomeCopy(exomecounts["chr1"],sample.name="HG00731.het.dup",
                   X.names=c("bg","GC","GC.sq","width"),S=0:6,d=2)
  show(fit)


###################################################
### code chunk number 13: exomeCopy.Rnw:242-243
###################################################
  copyCountSegments(fit)


###################################################
### code chunk number 14: exomeCopy.Rnw:249-252
###################################################
  cols <- c("red","orange","black","deepskyblue","blue","blue2","blue4")
  plot(fit,col=cols)
  abline(v=c(start(bounds),end(bounds)))


###################################################
### code chunk number 15: exomeCopy.Rnw:260-266
###################################################
runExomeCopy <- function(idx,rdata,seq.loop,sample.loop) {
  seq.name <- seq.loop[idx]
  sample.name <- sample.loop[idx]
  exomeCopy(rdata[seq.name],sample.name,
            X.names=c("bg","GC","GC.sq","width"),S=0:6,d=2)
}


###################################################
### code chunk number 16: exomeCopy.Rnw:271-279
###################################################
seqs <- c("chr1")
samples <- paste(exome.samples[1:4],cnv.type,sep=".")
nseqs <- length(seqs)
nsamples <- length(samples)
seq.loop <- rep(seqs,times=nsamples)
sample.loop <- rep(samples,each=nseqs)
fit.list <- lapply(seq_len(nseqs*nsamples),
                   runExomeCopy,exomecounts,seq.loop,sample.loop)


###################################################
### code chunk number 17: exomeCopy.Rnw:285-290
###################################################
par(mfrow=c(4,1),mar=c(3,3,1,1))
for (i in 1:4) {
  plot(fit.list[[i]],main="",xlab="",ylab="",show.legend=FALSE,col=cols)
  abline(v=c(start(bounds),end(bounds)))
}


###################################################
### code chunk number 18: exomeCopy.Rnw:298-300
###################################################
fit.list[[1]]@init.par$beta.hat
fit.list[[1]]@final.par$beta


###################################################
### code chunk number 19: exomeCopy.Rnw:306-313
###################################################
par(mfrow=c(2,2),mar=c(3,5,1,1))
for (i in 1:4) {
  barplot(rbind(fit.list[[i]]@final.par$beta,fit.list[[i]]@init.par$beta.hat),
          horiz=TRUE,las=1,beside=TRUE,col=c("white","darkgrey"),
          xlim=c(-150,350))
  legend("topright",legend=c("initial","final"),fill=c("darkgrey","white"))
}


###################################################
### code chunk number 20: exomeCopy.Rnw:319-320
###################################################
sessionInfo()


