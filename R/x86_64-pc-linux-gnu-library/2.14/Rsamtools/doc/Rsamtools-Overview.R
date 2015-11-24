### R code from vignette source 'Rsamtools-Overview.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(width=60)


###################################################
### code chunk number 2: preliminaries
###################################################
library(Rsamtools)


###################################################
### code chunk number 3: ScanBamParam
###################################################
which <- RangesList(seq1=IRanges(1000, 2000),
                    seq2=IRanges(c(100, 1000), c(1000, 2000)))
what <- c("rname", "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(which=which, what=what)


###################################################
### code chunk number 4: scanBam
###################################################
bamFile <- 
    system.file("extdata", "ex1.bam", package="Rsamtools")
bam <- scanBam(bamFile, param=param)


###################################################
### code chunk number 5: scanBam-elts-1
###################################################
class(bam)
names(bam)


###################################################
### code chunk number 6: scanBam-elts-2
###################################################
class(bam[[1]])
names(bam[[1]])


###################################################
### code chunk number 7: scanBam-elts-class
###################################################
sapply(bam[[1]], class)


###################################################
### code chunk number 8: scanBam-to-IRanges
###################################################
lst <- lapply(names(bam[[1]]), function(elt) {
    do.call(c, unname(lapply(bam, "[[", elt)))
})
names(lst) <- names(bam[[1]])


###################################################
### code chunk number 9: lst-to-DataFrame
###################################################
head(do.call("DataFrame", lst))


###################################################
### code chunk number 10: indexed-file
###################################################
list.files(dirname(bamFile), pattern="ex1.bam(.bai)?")


###################################################
### code chunk number 11: bam-remote (eval = FALSE)
###################################################
## which <- RangesList("6"=IRanges(100000L, 110000L))
## param <- ScanBamParam(which=which, what=scanBamWhat())
## na19240bam <- scanBam(na19240url, param=param)


###################################################
### code chunk number 12: summaryFunction
###################################################
summaryFunction <- 
    function(seqname, bamFile, ...)
{
    param <- ScanBamParam(what=c('pos', 'qwidth'),
                          which=GRanges(seqname, IRanges(1, 1e7)),
                          flag=scanBamFlag(isUnmappedQuery=FALSE))
    x <- scanBam(bamFile, ..., param=param)[[1]]
    coverage(IRanges(x[["pos"]], width=x[["qwidth"]]))
}


###################################################
### code chunk number 13: summaryByChromosome
###################################################
seqnames <- paste("seq", 1:2, sep="")
cvg <- lapply(seqnames, summaryFunction, bamFile)
names(cvg) <- seqnames
cvg


###################################################
### code chunk number 14: BamViews-parts
###################################################
library(GenomicFeatures)
bamRanges <- local({
    fl <- system.file("extdata", "CaffeineTxdb.sqlite", 
                      package="Rsamtools")
    transcripts(loadFeatures(fl))
})
slxMaq09 <- local({
    fl <- system.file("extdata", "slxMaq09_urls.txt", 
                      package="Rsamtools")
    readLines(fl)
})


###################################################
### code chunk number 15: BamViews-construct
###################################################
bamExperiment <- 
    list(description="Caffeine metabolism views on 1000 genomes samples",
         created=date())
bv <- BamViews(slxMaq09, bamRanges=bamRanges, 
               bamExperiment=bamExperiment)
metadata(bamSamples(bv)) <- 
    list(description="Solexa/MAQ samples, August 2009",
         created="Thu Mar 25 14:08:42 2010")


###################################################
### code chunk number 16: BamViews-query
###################################################
bamExperiment(bv)


###################################################
### code chunk number 17: bamIndicies (eval = FALSE)
###################################################
## bamIndexDir <- tempfile()
## indexFiles <- paste(bamPaths(bv), "bai", sep=".")
## dir.create(bamIndexDir)
## idxFiles <- mapply(download.file, indexFiles,
##                    file.path(bamIndexDir, basename(indexFiles)) ,
##                    MoreArgs=list(method="curl"))


###################################################
### code chunk number 18: readBamGappedAlignments (eval = FALSE)
###################################################
## library(multicore)
## olaps <- readBamGappedAlignments(bv)


###################################################
### code chunk number 19: olaps
###################################################
load(system.file("extdata", "olaps.Rda", package="Rsamtools"))
olaps
head(olaps[[1]])


###################################################
### code chunk number 20: read-lengths
###################################################
head(t(sapply(olaps, function(elt) range(qwidth(elt)))))


###################################################
### code chunk number 21: focal
###################################################
rng <- bamRanges(bv)[1]
strand(rng) <- "*"
olap1 <- endoapply(olaps, subsetByOverlaps, rng)
olap1
head(olap1[[24]])


###################################################
### code chunk number 22: olap-cvg
###################################################
minw <- min(sapply(olap1, function(elt) min(start(elt))))
maxw <- max(sapply(olap1, function(elt) max(end(elt))))
cvg <- endoapply(olap1, coverage,
                 shift=-start(ranges(bamRanges[1])),
                 width=width(ranges(bamRanges[1])))
cvg
cvg[[1]]


###################################################
### code chunk number 23: olap-cvg-as-m
###################################################
m <- matrix(unlist(lapply(cvg, lapply, as.vector)),
            ncol=length(cvg))
summary(rowSums(m))
summary(colSums(m))


###################################################
### code chunk number 24: sessionInfo
###################################################
packageDescription("Rsamtools")
sessionInfo()


###################################################
### code chunk number 25: caffeine-kegg
###################################################
library(KEGG.db)
kid <- revmap(KEGGPATHID2NAME)[["Caffeine metabolism"]]
egid <- KEGGPATHID2EXTID[[sprintf("hsa%s", kid)]]


###################################################
### code chunk number 26: caffeine-txdb (eval = FALSE)
###################################################
## library(biomaRt)
## mart <- useMart("ensembl", "hsapiens_gene_ensembl")
## ensid <- getBM(c("ensembl_transcript_id"), filters="entrezgene",
##                values=egid, mart=mart)[[1]]
## library(GenomicFeatures)
## txdb <- makeTranscriptDbFromBiomart(transcript_ids=ensid)


###################################################
### code chunk number 27: caffeine-txdb-load
###################################################
library(GenomicFeatures)
fl <- system.file("extdata", "CaffeineTxdb.sqlite", 
                  package="Rsamtools")
txdb <- loadFeatures(fl)


###################################################
### code chunk number 28: txdb-transcripts
###################################################
bamRanges <- transcripts(txdb)


###################################################
### code chunk number 29: bam-avail (eval = FALSE)
###################################################
## library(RCurl)
## ftpBase <- 
##     "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/data/"
## indivDirs <- 
##     strsplit(getURL(ftpBase, ftplistonly=TRUE), "\n")[[1]]
## alnDirs <- 
##     paste(ftpBase, indivDirs, "/alignment/", sep="")
## urls0 <- 
##     strsplit(getURL(alnDirs, dirlistonly=TRUE), "\n")


###################################################
### code chunk number 30: bam-index (eval = FALSE)
###################################################
## urls <- urls[sapply(urls0, length) != 0]
## fls0 <- unlist(unname(urls0))
## fls1 <- fls0[grepl("bai$", fls0)]
## fls <- fls1[sapply(strsplit(fls1, "\\."), length)==7]


###################################################
### code chunk number 31: slxMaq09 (eval = FALSE)
###################################################
## urls1 <- 
##     Filter(function(x) length(x) != 0,
##            lapply(urls, 
##                   function(x) x[grepl("SLX.maq.*2009_08.*bai$", x)]))
## slxMaq09.bai <- 
##    mapply(paste, names(urls1), urls1, sep="", USE.NAMES=FALSE)
## slxMaq09 <- sub(".bai$", "", slxMaq09.bai) #$


###################################################
### code chunk number 32: bamIndicies (eval = FALSE)
###################################################
## bamIndexDir <- tempfile()
## dir.create(bamIndexDir)
## idxFiles <- mapply(download.file, slxMaq09.bai, 
##                    file.path(bamIndexDir, basename(slxMaq09.bai)) ,
##                    MoreArgs=list(method="curl"))


