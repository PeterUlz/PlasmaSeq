### R code from vignette source 'BSgenomeForge.Rnw'

###################################################
### code chunk number 1: BSgenomeForge.Rnw:162-165
###################################################
library(Biostrings)
file <- system.file("extdata", "ce2chrM.fa", package="BSgenome")
fasta.info(file)


###################################################
### code chunk number 2: BSgenomeForge.Rnw:494-499
###################################################
library(BSgenome)
seed_files <- system.file("extdata", "GentlemanLab", package="BSgenome")
list.files(seed_files, pattern="-seed$")
rn4_seed <- list.files(seed_files, pattern="rn4", full.names=TRUE)
cat(readLines(rn4_seed), sep="\n")


###################################################
### code chunk number 3: BSgenomeForge.Rnw:512-514 (eval = FALSE)
###################################################
## library(BSgenome)
## forgeBSgenomeDataPkg("path/to/my/seed")


###################################################
### code chunk number 4: BSgenomeForge.Rnw:545-546
###################################################
sessionInfo()


