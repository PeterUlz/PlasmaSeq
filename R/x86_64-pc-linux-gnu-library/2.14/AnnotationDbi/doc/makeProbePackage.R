### R code from vignette source 'makeProbePackage.Rnw'

###################################################
### code chunk number 1: startup
###################################################
library("AnnotationDbi")


###################################################
### code chunk number 2: makeprobepackage
###################################################
filename <- system.file("extdata", "HG-U95Av2_probe_tab.gz", 
                          package="AnnotationDbi")
outdir   <- tempdir()
me       <- "Wolfgang Huber <w.huber@dkfz.de>"
species  <- "Homo_sapiens"
makeProbePackage("HG-U95Av2",
                 datafile   = gzfile(filename, open="r"),
                 outdir     = outdir,
                 maintainer = me,
                 species    = species,
                 version    = "0.0.1", 
                 force      = TRUE)
dir(outdir)


