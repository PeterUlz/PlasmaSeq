## A script for Marina Granovskaia to make Along-Chromosome Plots
## (C) W. Huber 2006
##
## Please specify the parameters in the following lines:

## Name of RNA and of DNA hybes.
## RNA must be a single file, for DNA multiple files can be separated by space
rna = "050208_mRNA30minx4_RH6.cel.gz"
## rna = "050209_mRNAx4_30min_re-hybe_RH6.cel.gz"
## rna = "050218_polyA-RNA_RH6_4x15min.cel.gz"

## rna = "05_04_27_2xpolyA_NAP3.cel.gz"

dna = "09_11_04_S96_genDNA_16hrs_45C_noDMSO.cel.gz  041119_S96genDNA_re-hybe.cel.gz   041120_S96genDNA_re-hybe.cel.gz"

## Set this to FALSE for newer scanner software, TRUE for older
## scanner software (that still rotated the .DAT files)
rotated = TRUE

## Directory Path to the CEL files
celpath =  "/ebi/research/huber/Projects/tilingArray/Celfiles"

## Which region should be plotted? (chromosome, start and end coordinates)
chr   = 10
start = 550000
end   = 580000
  
## ------------------------------------------------------------------------
## Now we are done with setting the parameters and can run the code
## ------------------------------------------------------------------------
library("tilingArray", verbose=FALSE)
library("davidTiling", verbose=FALSE)  ## Get these from Bioconductor

dna = strsplit(dna, " +")[[1]]
fnpretty = sub(".cel.*g*z*", "", rna, ignore.case=TRUE)

## Step 1: read the cel files
a = readCel2eSet(c(rna, dna), path=celpath, rotated=rotated)

## Step 2: DNA normalization
cat("Normalizing data")
data("probeAnno")
data("gff")

isRNA = a$filename %in% rna
isDNA = a$filename %in% dna
x = normalizeByReference(a[,isRNA], a[,isDNA], 
  pm=PMindex(probeAnno), background=BGindex(probeAnno),
  plotFileNames=sprintf("normalization-background-%s.pdf", fnpretty))

## Step 3: plot along chromosome
cat(".\nWriting the PDF file")
chr=as.integer(chr)
start=as.integer(start)
end=as.integer(end)
outfilename = sprintf("%s--chr%02d:%d-%d.pdf", fnpretty[1], chr, start, end)

pdf(file=outfilename, width=12, height=9)
trsf = function(x) { rg=range(x, na.rm=TRUE); log2(x-rg[1]+0.05*diff(rg)) }
plotAlongChrom(y=trsf(exprs(x)), probeAnno=probeAnno, gff=gff, chr=chr,
               coord=c(start, end), doLegend=TRUE, main=fnpretty[1])

dev.off()
cat(".\n\n")
