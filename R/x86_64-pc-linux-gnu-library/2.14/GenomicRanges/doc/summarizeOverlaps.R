### R code from vignette source 'summarizeOverlaps.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(width=72)


###################################################
### code chunk number 2: firstExample
###################################################
library(Rsamtools)
library(DESeq)
library(edgeR)

fls = list.files(system.file("extdata",package="GenomicRanges"),
    recursive=TRUE, pattern="*bam$", full=TRUE)
bfl <- BamFileList(fls)

features <- GRanges(
    seqnames = Rle(c("chr2L", "chr2R", "chr2L", "chr2R", "chr2L", "chr2R",
        "chr2L", "chr2R", "chr2R", "chr3L", "chr3L")),
    strand = strand(rep("+", 11)),
    ranges = IRanges(start=c(1000, 2000, 3000, 3600, 7000, 7500, 4000, 4000, 
        3000, 5000, 5400), width=c(500, 900, 500, 300, 600, 300, 500, 900, 500, 
        500, 500))
)

olap <- summarizeOverlaps(features, bfl)

deseq <- newCountDataSet(countData=assays(olap)$counts,  
                           conditions=rownames(colData(olap)))

edger <- DGEList(counts=assays(olap)$counts, group=rownames(colData(olap)))


###################################################
### code chunk number 3: simple
###################################################
rd <- GappedAlignments("a", rname = Rle("chr1"), pos = as.integer(100),
    cigar = "300M", strand = strand("+"))

gr1 <- GRanges("chr1", IRanges(start=50, width=150), strand="+")
gr2 <- GRanges("chr1", IRanges(start=350, width=150), strand="+")


###################################################
### code chunk number 4: simpleGRanges
###################################################
gr <- c(gr1, gr2)
data.frame(union = assays(summarizeOverlaps(gr, rd))$counts,
           intStrict = assays(summarizeOverlaps(gr, rd,
               mode="IntersectionStrict"))$counts,
           intNotEmpty = assays(summarizeOverlaps(gr, rd,
               mode="IntersectionNotEmpty"))$counts)


###################################################
### code chunk number 5: simpleGRangesList
###################################################
grl <- GRangesList(c(gr1, gr2))
data.frame(union = assays(summarizeOverlaps(grl, rd))$counts,
           intStrict = assays(summarizeOverlaps(grl, rd,
               mode="IntersectionStrict"))$counts,
           intNotEmpty = assays(summarizeOverlaps(grl, rd,
               mode="IntersectionNotEmpty"))$counts)


###################################################
### code chunk number 6: data
###################################################
group_id <- c("A", "B", "C", "C", "D", "D", "E", "F", "G", "G", "H", "H")
features <- GRanges(
    seqnames = Rle(c("chr1", "chr2", "chr1", "chr1", "chr2", "chr2",
        "chr1", "chr1", "chr2", "chr2", "chr1", "chr1")),
    strand = strand(rep("+", length(group_id))),
    ranges = IRanges(
        start=c(1000, 2000, 3000, 3600, 7000, 7500, 4000, 4000, 3000, 3350, 5000, 5400),
        width=c(500, 900, 500, 300, 600, 300, 500, 900, 150, 200, 500, 500)),
   DataFrame(group_id)
)

reads <- GappedAlignments(
    names = c("a","b","c","d","e","f","g"),
    rname = Rle(c(rep(c("chr1", "chr2"), 3), "chr1")),
    pos = as.integer(c(1400, 2700, 3400, 7100, 4000, 3100, 5200)),
    cigar = c("500M", "100M", "300M", "500M", "300M", "50M200N50M", "50M150N50M"),
    strand = strand(rep.int("+", 7L)))



###################################################
### code chunk number 7: GRanges
###################################################
data.frame(union = assays(summarizeOverlaps(features, reads))$counts,
           intStrict = assays(summarizeOverlaps(features, reads,
               mode="IntersectionStrict"))$counts,
           intNotEmpty = assays(summarizeOverlaps(features, reads,
               mode="IntersectionNotEmpty"))$counts)


###################################################
### code chunk number 8: lst
###################################################
lst <- split(features, values(features)[["group_id"]])
length(lst)


###################################################
### code chunk number 9: GRangesList
###################################################
data.frame(union = assays(summarizeOverlaps(lst, reads))$counts,
           intStrict = assays(summarizeOverlaps(lst, reads,
               mode="IntersectionStrict"))$counts,
           intNotEmpty = assays(summarizeOverlaps(lst, reads,
               mode="IntersectionNotEmpty"))$counts)


###################################################
### code chunk number 10: pasilla_features
###################################################
library(pasilla)
library(rtracklayer)
library(Rsamtools)

gff <- import(system.file("extdata", "Dmel.BDGP5.25.62.DEXSeq.chr.gff",
    package = "pasilla"), "gff1")
features <- as(gff, "GRanges")
head(features[,1])


###################################################
### code chunk number 11: pasilla_exons
###################################################
exons <- features[values(features)[["type"]] == "exonic_part"]

st <- strsplit(gsub("\"", "", values(exons)[["group"]]), ";") 
exonID <- do.call(c,
              lapply(st, function(x) {
                 gsub("[^0-9]", "", x[2])}))

geneID <- do.call(c,
              lapply(st, function(x) {
                 gsub(" gene_id ", "", x[3])}))


###################################################
### code chunk number 12: pasilla_param
###################################################
param <- ScanBamParam(
             what='qual',
             which=GRanges("chr2L", IRanges(1, 1e+6)),
             flag=scanBamFlag(isUnmappedQuery=FALSE, isPaired=NA))
bamTag(param) <- "NH" 


###################################################
### code chunk number 13: pasilla_count (eval = FALSE)
###################################################
## fls <- c("treated1.bam", "untreated1.bam", "untreated2.bam") 
## path <- "pathToBAMFiles"
## bamFiles <- BamFileList(file.path(paste(path, fls, sep=""))) 
## se_exons <- summarizeOverlaps(exons, bamFiles, mode="Union")


###################################################
### code chunk number 14: pasilla_exoncountset (eval = FALSE)
###################################################
## library(DEXSeq)
## expdata = new("MIAME",
##               name="pasilla knockdown",
##               lab="Genetics and Developmental Biology, University of 
##                   Connecticut Health Center",
##               contact="Dr. Brenton Graveley",
##               title="modENCODE Drosophila pasilla RNA Binding Protein RNAi 
##                   knockdown RNA-Seq Studies",
##               url="http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE18508",
##               abstract="RNA-seq of 3 biological replicates of from the Drosophila
##                   melanogaster S2-DRSC cells that have been RNAi depleted of mRNAs 
##                   encoding pasilla, a mRNA binding protein and 4 biological replicates 
##                   of the the untreated cell line.")
##               pubMedIds(expdata) <- "20921232"
## 
## design <- data.frame(
##               condition=c("treated", "untreated", "untreated"),
##               replicate=c(1,1,2),
##               type=rep("single-read", 3),
##               countfiles=colData(se_exons)[,1], stringsAsFactors=TRUE)
## 
## pasillaECS <- newExonCountSet(
##                   countData=assays(se_exons)$counts,
##                   design=design,
##                   exonIDs=factor(exonID), 
##                   geneIDs=factor(geneID))
## 
## experimentData(pasillaECS) <- expdata
## sampleNames(pasillaECS) = colnames(se_exons)


###################################################
### code chunk number 15: pasilla_bingenes (eval = FALSE)
###################################################
## genetable = geneCountTable(pasillaECS)
## pasillaCDS = newCountDataSet(countData=genetable, conditions=design)
## experimentData(pasillaCDS) = expdata


###################################################
### code chunk number 16: pasilla_genes (eval = FALSE)
###################################################
## genes <- features[values(features)[["type"]] == "aggregate_gene"]
## se_genes <- summarizeOverlaps(genes, bamFiles, mode="Union")
## 
## pasillaCDS_alt <- newCountDataSet(countData=assays(se_genes)$counts,
##     conditions=design) 
## experimentData(pasillaCDS_alt) = expdata


