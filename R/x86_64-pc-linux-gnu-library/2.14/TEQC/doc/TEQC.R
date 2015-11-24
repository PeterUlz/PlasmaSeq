### R code from vignette source 'TEQC.Rnw'

###################################################
### code chunk number 1: load
###################################################
library(TEQC)

exptPath <- system.file("extdata", package="TEQC")
targets <- get.targets(targetsfile=paste(exptPath, "ExampleSet_Targets.bed", sep="/"), chrcol=1, startcol=2, endcol=3, skip=0)
targets


###################################################
### code chunk number 2: fraction.target
###################################################
ft <- fraction.target(targets, genome="hg19")
ft


###################################################
### code chunk number 3: reads
###################################################
reads <- get.reads(paste(exptPath, "ExampleSet_Reads.bed", sep="/"), chrcol=1, startcol=2, endcol=3, idcol=4, zerobased=F, skip=0)
reads


###################################################
### code chunk number 4: reads2pairs
###################################################
readpairs <- reads2pairs(reads)
readpairs


###################################################
### code chunk number 5: insertsizehist
###################################################
insert.size.hist(readpairs, breaks=10)


###################################################
### code chunk number 6: chrombarplot
###################################################
chrom.barplot(reads, targets)


###################################################
### code chunk number 7: fraction.reads.target
###################################################
fr <- fraction.reads.target(reads, targets)
fr


###################################################
### code chunk number 8: withoffset
###################################################
fraction.reads.target(reads, targets, Offset=100)


###################################################
### code chunk number 9: enrichment
###################################################
fr / ft


###################################################
### code chunk number 10: fracpairs
###################################################
fraction.reads.target(readpairs, targets)


###################################################
### code chunk number 11: coverage.target
###################################################
Coverage <- coverage.target(reads, targets, perTarget=T, perBase=T)
Coverage
targets2 <- Coverage$targetCoverages


###################################################
### code chunk number 12: readspertarget
###################################################
targets2 <- readsPerTarget(reads, targets2)
targets2


###################################################
### code chunk number 13: covered.k
###################################################
covered.k(Coverage$coverageTarget, k=c(1, 5, 10))


###################################################
### code chunk number 14: coveragehist
###################################################
coverage.hist(Coverage$coverageTarget, covthreshold=8)


###################################################
### code chunk number 15: coverageuniformity
###################################################
coverage.uniformity(Coverage)


###################################################
### code chunk number 16: coveragetargetlength
###################################################
par(mfrow=c(1,2))
coverage.targetlength.plot(targets2, plotcolumn="nReads", pch=16, cex=1.5)
coverage.targetlength.plot(targets2, plotcolumn="avgCoverage", pch=16, cex=1.5)


###################################################
### code chunk number 17: get.baits
###################################################
baitsfile <- paste(exptPath, "ExampleSet_Baits.txt", sep="/")
baits <- get.baits(baitsfile, chrcol=3, startcol=4, endcol=5, seqcol=2)


###################################################
### code chunk number 18: coverageGC
###################################################
coverage.GC(Coverage$coverageAll, baits, pch=16, cex=1.5)


###################################################
### code chunk number 19: coverageplot
###################################################
coverage.plot(Coverage$coverageAll, targets, Offset=100, chr="chr1", Start=11157524, End=11158764)


###################################################
### code chunk number 20: make.wigfiles (eval = FALSE)
###################################################
## make.wigfiles(Coverage$coverageAll)


###################################################
### code chunk number 21: duplicatesbarplot
###################################################
duplicates.barplot(reads, targets)


###################################################
### code chunk number 22: duplicatesbarplot2
###################################################
duplicates.barplot(readpairs, targets, ylab="Fraction of read pairs")


###################################################
### code chunk number 23: collapse1
###################################################
dupfun <- function(x) duplicated(x$ranges)
params <- RDApplyParams(rangedData=reads, applyFun=dupfun)
dups <- unlist(rdapply(params))
reads.collapsed <- reads[!dups,,drop=T]


###################################################
### code chunk number 24: collapse2
###################################################
params2 <- RDApplyParams(rangedData=readpairs, applyFun=dupfun)
dups2 <- unlist(rdapply(params2))
ID.nondups <- readpairs$ID[!dups2]
sel <- reads$ID %in% ID.nondups
reads.collapsed.pairs <- reads[sel,,drop=T]


###################################################
### code chunk number 25: cov1
###################################################
coverage.target(reads.collapsed, targets, perBase=F, perTarget=F)


###################################################
### code chunk number 26: cov2
###################################################
coverage.target(reads.collapsed.pairs, targets, perBase=F, perTarget=F)


###################################################
### code chunk number 27: newsample
###################################################
r <- sample(nrow(reads), 0.1 * nrow(reads))
reads2 <- reads[-r,,drop=T]
Coverage2 <- coverage.target(reads2, targets, perBase=T)


###################################################
### code chunk number 28: coveragedensity
###################################################
covlist <- list(Coverage, Coverage2)
par(mfrow=c(1,2))
coverage.density(covlist)
coverage.density(covlist, normalized=F)


###################################################
### code chunk number 29: coverageuniformity2
###################################################
coverage.uniformity(Coverage, addlines=F)
coverage.uniformity(Coverage2, addlines=F, add=T, col="blue", lty=2)


###################################################
### code chunk number 30: coverageplot2
###################################################
coverage.plot(Coverage$coverageAll, targets, Offset=100, chr="chr1", Start=11157524, End=11158764)
coverage.plot(Coverage2$coverageAll, add=T, col.line=2, chr="chr1", Start=11157524, End=11158764)


###################################################
### code chunk number 31: covcor
###################################################
coverage.correlation(covlist, plotfrac=0.1, cex.pch=4)


###################################################
### code chunk number 32: sessioninfo
###################################################
toLatex(sessionInfo())


