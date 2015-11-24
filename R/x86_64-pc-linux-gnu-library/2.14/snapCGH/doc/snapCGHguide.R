### R code from vignette source 'snapCGHguide.Rnw'

###################################################
### code chunk number 1: 1
###################################################
library(snapCGH)


###################################################
### code chunk number 2: 2
###################################################
datadir <- system.file("testdata", package="snapCGH")
targets <- readTargets("targets.txt", path=datadir)
RG1 <- read.maimages(targets$FileName, path=datadir, source = "genepix")


###################################################
### code chunk number 3: 3
###################################################
RG1 <- read.clonesinfo("cloneinfo.txt", RG1, path=datadir)
RG1$printer <- getLayout(RG1$genes)
types <- readSpotTypes("SpotTypes.txt", path=datadir)
RG1$genes$Status <- controlStatus(types, RG1)



###################################################
### code chunk number 4: 3a
###################################################
RG1$design <- c(-1,-1)


###################################################
### code chunk number 5: 4
###################################################
RG2 <- backgroundCorrect(RG1, method="minimum")


###################################################
### code chunk number 6: 5
###################################################
MA <- normalizeWithinArrays(RG2, method="median")


###################################################
### code chunk number 7: 6
###################################################
MA2 <- processCGH(MA,method.of.averaging=mean, ID = "ID")


###################################################
### code chunk number 8: segmentation1
###################################################
SegInfo.Hom <- runHomHMM(MA2, criteria = "AIC")


###################################################
### code chunk number 9: segmentation2 (eval = FALSE)
###################################################
## SegInfo.GLAD <- runGLAD(MA2)
## SegInfo.DNAcopy <- runDNAcopy(MA2)
## SegInfo.TilingArray <- runTilingArray(MA2)


###################################################
### code chunk number 10: segmentation3 (eval = FALSE)
###################################################
## SegInfo.Bio <- runBioHMM(MA2)


###################################################
### code chunk number 11: segmentation4
###################################################
SegInfo.Hom.merged <- mergeStates(SegInfo.Hom, MergeType = 1)


###################################################
### code chunk number 12: plotting1
###################################################
genomePlot(MA2, array = 1)                                               


###################################################
### code chunk number 13: plotting2
###################################################
genomePlot(MA2, array = 1)                                               


###################################################
### code chunk number 14: plotting3
###################################################
genomePlot(MA2, array = 1, chrom.to.plot = 8)


###################################################
### code chunk number 15: plotting4
###################################################
genomePlot(MA2, array = 1, chrom.to.plot = 8)


###################################################
### code chunk number 16: 14
###################################################
plotSegmentedGenome(SegInfo.Hom.merged, array = 1)


###################################################
### code chunk number 17: 14a
###################################################
plotSegmentedGenome(SegInfo.Hom.merged, array = 1)


###################################################
### code chunk number 18: 15
###################################################
Seg.DNAcopy <- runDNAcopy(MA2)
SegInfo.DNAcopy.merged <- mergeStates(Seg.DNAcopy)
plotSegmentedGenome(SegInfo.DNAcopy.merged, SegInfo.Hom.merged, array = 1,
                    chrom.to.plot = 1, colors = c("blue", "green"))


###################################################
### code chunk number 19: 15a
###################################################
Seg.DNAcopy <- runDNAcopy(MA2)
SegInfo.DNAcopy.merged <- mergeStates(Seg.DNAcopy)
plotSegmentedGenome(SegInfo.DNAcopy.merged, SegInfo.Hom.merged, array = 1,
                    chrom.to.plot = 1, colors = c("blue", "green"))


###################################################
### code chunk number 20: iplotting1 (eval = FALSE)
###################################################
## zoomGenome(SegInfo.Hom.merged, array = 1)


###################################################
### code chunk number 21: iplotting2 (eval = FALSE)
###################################################
## zoomChromosome(SegInfo.Hom.merged, array = 1, chrom.to.plot = 8)


###################################################
### code chunk number 22: Simulation1
###################################################
simulation <- simulateData(nArrays = 4)


###################################################
### code chunk number 23: Simulation2
###################################################
Sim.HomHMM <- runHomHMM(simulation)
Sim.DNAcopy <- runDNAcopy(simulation)
rates <- compareSegmentations(simulation, offset = 0, Sim.HomHMM, Sim.DNAcopy)


###################################################
### code chunk number 24: Simulation3
###################################################
rates


###################################################
### code chunk number 25: Simulation4
###################################################
par(mfrow = c(1,2))
boxplot(rates$TPR ~ row(rates$TPR), col = c("red", "blue"), main = "True Positive Rate")
boxplot(rates$FDR ~ row(rates$FDR), col = c("red", "blue"), main = "False Discovery Rate")


