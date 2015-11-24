# R script to create a tab-delimited flat file containing the probe-wise
## normalized intensities:

## 0. load libraries and data
library(tilingArray)
library(davidTiling)
data(davidTiling)
data(probeAnno)

## where to write the flat file to?
outDir <- "~/temp"

## 1. determine perfect match (PM) probes and background (BG) probes
isPM = logical(nrow(exprs(davidTiling)))
for(j in probeAnno$probeReverse)
  isPM[ as.character(j) != ""] = TRUE
isBG = (probeAnno$probeReverse$no_feature=="no" &
        probeAnno$probeDirect$no_feature=="no")

sum(isPM)
(tab=table(isPM, isBG))
# note: isBG is a subset of isPM

## 2. determine sample type
isRNA = davidTiling$nucleicAcid %in% c("poly(A) RNA","total RNA")
isDNA = davidTiling$nucleicAcid %in% "genomic DNA"
stopifnot(sum(isRNA)==5, sum(isDNA)==3)

## 3. do normalization:
xn2 = normalizeByReference(davidTiling[,isRNA], davidTiling[,isDNA], 
  pm=isPM, background=isBG)
#save(xn2, file="~/temp/xn2.rda")

## 4. write data frame:
arePolyA <- grep("poly",sampleNames(xn2))
areTotal <- grep("totcDNA",sampleNames(xn2))

resultDf <- data.frame()

for (chr in 1:17) { # for all chromosomes
  for (strand in c("+","-")) { # and both strands:
    cat(chr,strand,"...\n")
    # get index of probes on chr and strand:
    thisIdx <- get(paste(chr,strand,"index",sep="."),env=probeAnno)
    nTheseProbes <- length(thisIdx)
    thisPol <- rowMeans(exprs(xn2)[thisIdx,arePolyA])
    thisTot <- rowMeans(exprs(xn2)[thisIdx,areTotal])
    thischrstrDf <- data.frame(Chromosome=I(rep(chr, nTheseProbes)),
                               Strand=I(rep(strand, nTheseProbes)),
                               Start=I(get(paste(chr,strand,"start",sep="."),env=probeAnno)),
                               End=I(get(paste(chr,strand,"end",sep="."),env=probeAnno)),
                               Unique=I(get(paste(chr,strand,"unique",sep="."),env=probeAnno)),
                               PolyA=I(thisPol), Total=I(thisTot))
    resultDf <- rbind(resultDf, thischrstrDf)
  } # for both strands
}# for all chromosomes

# only report those genomic locations uniquely hit by a single probe:
areUnique <- (resultDf$Unique==0)
resultDf <- resultDf[areUnique,]

probeMiddle <- round((resultDf$Start + resultDf$End)/2)
resultDf$Position <- probeMiddle

probeOrder <- order(resultDf$Chromosome, resultDf$Strand, resultDf$Position)

resultDf2 <- resultDf[probeOrder,c("Chromosome","Strand","Position","PolyA","Total")]
resultDf2$PolyA <- round(resultDf2$PolyA, digits=4)
resultDf2$Total <- round(resultDf2$Total, digits=4)

write.table(resultDf2, file=file.path(outDir,"normProbeIntensities.txt"), sep="\t", quote=FALSE,
            row.names=FALSE, col.names=TRUE)

setwd(outDir)
system("zip -r9 normProbeIntensities.zip normProbeIntensities.txt")
