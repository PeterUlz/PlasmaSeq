### R code from vignette source 'GenomeSearching.Rnw'

###################################################
### code chunk number 1: b1
###################################################
library(BSgenome.Celegans.UCSC.ce2)
ls("package:BSgenome.Celegans.UCSC.ce2")
Celegans


###################################################
### code chunk number 2: b2
###################################################
class(Celegans)


###################################################
### code chunk number 3: b3
###################################################
organism(Celegans)
provider(Celegans)
providerVersion(Celegans)
seqnames(Celegans)
mseqnames(Celegans)


###################################################
### code chunk number 4: b4
###################################################
Celegans$chrI


###################################################
### code chunk number 5: b5
###################################################
chrI <- Celegans$chrI
length(chrI)


###################################################
### code chunk number 6: b6
###################################################
afI <- alphabetFrequency(chrI)
afI
sum(afI) == length(chrI)


###################################################
### code chunk number 7: b7
###################################################
p1 <- "ACCCAGGGC"
countPattern(p1, chrI)


###################################################
### code chunk number 8: b8
###################################################
countPattern(p1, chrI, max.mismatch=1)


###################################################
### code chunk number 9: b9
###################################################
m1 <- matchPattern(p1, chrI, max.mismatch=1)
m1[4:6]
class(m1)


###################################################
### code chunk number 10: b10
###################################################
mismatch(p1, m1[4:6])


###################################################
### code chunk number 11: b11
###################################################
p2 <- DNAString("AAGCCTAAGCCTAAGCCTAA")
m2 <- matchPattern(p2, chrI, max.mismatch=2)
m2[1:4]
p2 == m2[1:4]
mismatch(p2, m2[1:4])


###################################################
### code chunk number 12: b12
###################################################
m2[p2 == m2]
m2[p2 != m2]


###################################################
### code chunk number 13: c1
###################################################
ce2dict0_file <- system.file("extdata", "ce2dict0.fa", package="BSgenome")
ce2dict0 <- read.DNAStringSet(ce2dict0_file, "fasta")
ce2dict0


###################################################
### code chunk number 14: c2
###################################################
writeHits <- function(seqname, matches, strand, file="", append=FALSE)
{
    if (file.exists(file) && !append)
        warning("existing file ", file, " will be overwritten with 'append=FALSE'")
    if (!file.exists(file) && append)
        warning("new file ", file, " will have no header with 'append=TRUE'")
    hits <- data.frame(seqname=rep.int(seqname, length(matches)),
                       start=start(matches),
                       end=end(matches),
                       strand=rep.int(strand, length(matches)),
                       patternID=names(matches),
                       check.names=FALSE)
    write.table(hits, file=file, append=append, quote=FALSE, sep="\t",
                row.names=FALSE, col.names=!append)
}

runAnalysis1 <- function(dict0, outfile="")
{
    library(BSgenome.Celegans.UCSC.ce2)
    seqnames <- seqnames(Celegans)
    seqnames_in1string <- paste(seqnames, collapse=", ")
    cat("Target:", providerVersion(Celegans),
        "chromosomes", seqnames_in1string, "\n")
    append <- FALSE
    for (seqname in seqnames) {
        subject <- Celegans[[seqname]]
        cat(">>> Finding all hits in chromosome", seqname, "...\n")
        for (i in seq_len(length(dict0))) {
            patternID <- names(dict0)[i]
            pattern <- dict0[[i]]
            plus_matches <- matchPattern(pattern, subject)
            names(plus_matches) <- rep.int(patternID, length(plus_matches))
            writeHits(seqname, plus_matches, "+", file=outfile, append=append)
            append <- TRUE
            rcpattern <- reverseComplement(pattern)
            minus_matches <- matchPattern(rcpattern, subject)
            names(minus_matches) <- rep.int(patternID, length(minus_matches))
            writeHits(seqname, minus_matches, "-", file=outfile, append=append)
        }
        cat(">>> DONE\n")
    }
}


###################################################
### code chunk number 15: c3
###################################################
runAnalysis1(ce2dict0, outfile="ce2dict0_ana1.txt")


###################################################
### code chunk number 16: c4
###################################################
hits1 <- read.table("ce2dict0_ana1.txt", header=TRUE)
nrow(hits1)


###################################################
### code chunk number 17: c5
###################################################
table(hits1$seqname)


###################################################
### code chunk number 18: c6
###################################################
hits1_table <- table(hits1$patternID)
hits1_table


###################################################
### code chunk number 19: c7
###################################################
hits1_table[hits1_table == max(hits1_table)] # pattern(s) with more hits


###################################################
### code chunk number 20: c8
###################################################
setdiff(names(ce2dict0), hits1$patternID) # pattern(s) with no hits


###################################################
### code chunk number 21: c9
###################################################
plotGenomeHits <- function(bsgenome, seqnames, hits)
{
    chrlengths <- seqlengths(bsgenome)[seqnames]
    XMAX <- max(chrlengths)
    YMAX <- length(seqnames)
    plot.new()
    plot.window(c(1, XMAX), c(0, YMAX))
    axis(1)
    axis(2, at=seq_len(length(seqnames)), labels=rev(seqnames), tick=FALSE, las=1)
    ## Plot the chromosomes
    for (i in seq_len(length(seqnames)))
        lines(c(1, chrlengths[i]), c(YMAX + 1 - i, YMAX + 1 - i), type="l")
    ## Plot the hits
    for (i in seq_len(nrow(hits))) {
        seqname <- hits$seqname[i]
        y0 <- YMAX + 1 - match(seqname, seqnames)
        if (hits$strand[i] == "+") {
            y <- y0 + 0.05
            col <- "red"
        } else {
            y <- y0 - 0.05
            col <- "blue"
        }
        lines(c(hits$start[i], hits$end[i]), c(y, y), type="l", col=col, lwd=3)
    }
}


###################################################
### code chunk number 22: c10 (eval = FALSE)
###################################################
## plotGenomeHits(Celegans, seqnames(Celegans), hits1)


###################################################
### code chunk number 23: d1
###################################################
library(hgu95av2probe)
tmpseq <- DNAStringSet(hgu95av2probe$sequence)
someStats <- function(v)
{
    GC <- DNAString("GC")
    CG <- DNAString("CG")
    sapply(seq_len(length(v)),
           function(i) {
               y <- v[[i]]
               c(alphabetFrequency(y)[1:4],
                 GC=countPattern(GC, y),
                 CG=countPattern(CG, y))
           }
    )
}
someStats(tmpseq[1:10])


###################################################
### code chunk number 24: f1
###################################################
library(BSgenome.Hsapiens.UCSC.hg19)
chrY <- Hsapiens$chrY
chrY
chrM <- Hsapiens$chrM
chrM


###################################################
### code chunk number 25: f2
###################################################
active(masks(chrY))["RM"] <- TRUE
chrY


###################################################
### code chunk number 26: f3
###################################################
active(masks(chrY)) <- FALSE
active(masks(chrY))["AGAPS"] <- TRUE
alphabetFrequency(unmasked(chrY))
alphabetFrequency(chrY)


###################################################
### code chunk number 27: f4
###################################################
as(chrY, "XStringViews")


###################################################
### code chunk number 28: f5
###################################################
gaps(as(chrY, "XStringViews"))


###################################################
### code chunk number 29: f6
###################################################
width(gaps(as(chrY, "XStringViews")))


###################################################
### code chunk number 30: f7
###################################################
gaps(chrY)
alphabetFrequency(gaps(chrY))


###################################################
### code chunk number 31: f8
###################################################
af0 <- alphabetFrequency(unmasked(chrY))
af1 <- alphabetFrequency(chrY)
af2 <- alphabetFrequency(gaps(chrY))
all(af0 == af1 + af2)


###################################################
### code chunk number 32: f9
###################################################
active(masks(chrY)) <- TRUE
af1 <- alphabetFrequency(chrY)
af1
gaps(chrY)
af2 <- alphabetFrequency(gaps(chrY))
af2
all(af0 == af1 + af2)


###################################################
### code chunk number 33: f10
###################################################
Ebox <- "CANNTG"
active(masks(chrY)) <- FALSE
countPattern(Ebox, chrY, fixed=FALSE)


###################################################
### code chunk number 34: f11
###################################################
countPattern(Ebox, chrY, fixed=c(pattern=FALSE,subject=TRUE))


###################################################
### code chunk number 35: f12
###################################################
active(masks(chrY))[c("AGAPS", "AMB")] <- TRUE
alphabetFrequency(chrY, baseOnly=TRUE)  # no ambiguities
countPattern(Ebox, chrY, fixed=FALSE)


###################################################
### code chunk number 36: f13
###################################################
chr2 <- Hsapiens$chr2
active(masks(chr2))[-2] <- FALSE
alphabetFrequency(gaps(chr2))


###################################################
### code chunk number 37: e1
###################################################
runOneStrandAnalysis <- function(dict0, bsgenome, seqnames, strand,
                                 outfile="", append=FALSE)
{
    cat("\nTarget: strand", strand, "of", providerVersion(bsgenome),
        "chromosomes", paste(seqnames, collapse=", "), "\n")
    if (strand == "-")
        dict0 <- reverseComplement(dict0)
    pdict <- PDict(dict0)
    for (seqname in seqnames) {
        subject <- bsgenome[[seqname]]
        cat(">>> Finding all hits in strand", strand, "of chromosome", seqname, "...\n")
        mindex <- matchPDict(pdict, subject)
        matches <- extractAllMatches(subject, mindex)
        writeHits(seqname, matches, strand, file=outfile, append=append)
        append <- TRUE
        cat(">>> DONE\n")
    }
}

runAnalysis2 <- function(dict0, outfile="")
{
    library(BSgenome.Celegans.UCSC.ce2)
    seqnames <- seqnames(Celegans)
    runOneStrandAnalysis(dict0, Celegans, seqnames, "+", outfile=outfile, append=FALSE)
    runOneStrandAnalysis(dict0, Celegans, seqnames, "-", outfile=outfile, append=TRUE)
}


###################################################
### code chunk number 38: e2
###################################################
ce2dict0cw15 <- DNAStringSet(ce2dict0, end=15)


###################################################
### code chunk number 39: e3
###################################################
runAnalysis2(ce2dict0cw15, outfile="ce2dict0cw15_ana2.txt")


###################################################
### code chunk number 40: g1
###################################################
sessionInfo()


