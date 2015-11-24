### R code from vignette source 'cn.mops.Rnw'

###################################################
### code chunk number 1: cn.mops.Rnw:40-44
###################################################
options(width=75)
set.seed(0)
library(cn.mops)
cn.mopsVersion <- packageDescription("cn.mops")$Version


###################################################
### code chunk number 2: cn.mops.Rnw:137-138
###################################################
library(cn.mops)


###################################################
### code chunk number 3: cn.mops.Rnw:147-149 (eval = FALSE)
###################################################
## BAMFiles <- list.files(pattern=".bam$")
## bamDataRanges <- getReadCountsFromBAM(BAMFiles)


###################################################
### code chunk number 4: cn.mops.Rnw:153-154 (eval = FALSE)
###################################################
## res <- cn.mops(bamDataRanges)


###################################################
### code chunk number 5: cn.mops.Rnw:159-160 (eval = FALSE)
###################################################
## plot(res,which=1)


###################################################
### code chunk number 6: cn.mops.Rnw:164-169
###################################################
data(cn.mops)
resCNMOPS <- cn.mops(XRanges)
pdf("003.pdf")
plot(resCNMOPS,which=7,toFile=TRUE)
dev.off()


###################################################
### code chunk number 7: cn.mops.Rnw:228-232
###################################################
BAMFiles <- list.files(system.file("extdata", package="cn.mops"),pattern=".bam$",
		full.names=TRUE)
bamDataRanges <- getReadCountsFromBAM(BAMFiles,
		sampleNames=paste("Sample",1:3))


###################################################
### code chunk number 8: cn.mops.Rnw:237-238
###################################################
(bamDataRanges)


###################################################
### code chunk number 9: cn.mops.Rnw:257-259
###################################################
data(cn.mops)
ls()


###################################################
### code chunk number 10: cn.mops.Rnw:264-265
###################################################
head(XRanges[,1:3])


###################################################
### code chunk number 11: cn.mops.Rnw:268-269 (eval = FALSE)
###################################################
## resCNMOPS <- cn.mops(XRanges)


###################################################
### code chunk number 12: cn.mops.Rnw:277-278
###################################################
head(X[,1:3])


###################################################
### code chunk number 13: cn.mops.Rnw:281-282 (eval = FALSE)
###################################################
## resCNMOPSX <- cn.mops(X)


###################################################
### code chunk number 14: cn.mops.Rnw:288-289 (eval = FALSE)
###################################################
## all(individualCall(resCNMOPSX)==individualCall(resCNMOPS))


###################################################
### code chunk number 15: cn.mops.Rnw:295-296 (eval = FALSE)
###################################################
## (resCNMOPS)


###################################################
### code chunk number 16: cn.mops.Rnw:300-301
###################################################
cnvs(resCNMOPS)[1:5]


###################################################
### code chunk number 17: cn.mops.Rnw:306-307
###################################################
cnvr(resCNMOPS)[1,1:5]


###################################################
### code chunk number 18: cn.mops.Rnw:313-314
###################################################
(CNVRanges[15,1:5])


###################################################
### code chunk number 19: cn.mops.Rnw:320-322
###################################################
ranges(cnvr(resCNMOPS))[1:2]
ranges(cnvr(resCNMOPS)) %in% ranges(CNVRanges)


###################################################
### code chunk number 20: cn.mops.Rnw:328-329 (eval = FALSE)
###################################################
## help(CNVDetectionResult)


###################################################
### code chunk number 21: cn.mops.Rnw:338-339 (eval = FALSE)
###################################################
## segplot(resCNMOPS,sampleIdx=13)


###################################################
### code chunk number 22: cn.mops.Rnw:342-345
###################################################
pdf("002.pdf")
segplot(resCNMOPS,sampleIdx=13)
dev.off()


###################################################
### code chunk number 23: cn.mops.Rnw:364-365 (eval = FALSE)
###################################################
## plot(resCNMOPS,which=1)


###################################################
### code chunk number 24: cn.mops.Rnw:368-371
###################################################
pdf("001.pdf")
plot(resCNMOPS,which=1,toFile=TRUE)
dev.off()


###################################################
### code chunk number 25: cn.mops.Rnw:434-439 (eval = FALSE)
###################################################
## BAMFiles <- list.files(pattern=".bam$")
## segments <- read.table("targetRegions.bed",sep="\t",as.is=TRUE)
## gr <- GRanges(segments[,1],IRanges(segments[,2],segments[,3]))
## X <- getSegmentReadCountsFromBAM(BAMFiles,GR=gr)
## resCNMOPS <- cn.mops(X)


###################################################
### code chunk number 26: cn.mops.Rnw:453-455
###################################################
XchrX <- normalizeChromosomes(X[1:500, ],ploidy=c(rep(1,10),rep(2,30)))
cnvr(cn.mops(XchrX,norm=FALSE))


###################################################
### code chunk number 27: cn.mops.Rnw:474-475 (eval = FALSE)
###################################################
## toBibtex(citation("cn.mops"))


