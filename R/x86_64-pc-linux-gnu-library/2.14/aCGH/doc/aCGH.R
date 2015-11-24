### R code from vignette source 'aCGH.Rnw'

###################################################
### code chunk number 1: aCGH.Rnw:58-73
###################################################
library(aCGH)

datadir <- system.file(package = "aCGH")
datadir <- paste(datadir, "/examples", sep="")

clones.info <-
      read.table(file = file.path(datadir, "clones.info.ex.txt"),
                 header = T, sep = "\t", quote="", comment.char="")
log2.ratios <-
      read.table(file = file.path(datadir, "log2.ratios.ex.txt"),
                 header = T, sep = "\t", quote="", comment.char="")
pheno.type <-
      read.table(file = file.path(datadir, "pheno.type.ex.txt"),
                 header = T, sep = "\t", quote="", comment.char="")
ex.acgh <- create.aCGH(log2.ratios, clones.info, pheno.type)


###################################################
### code chunk number 2: aCGH.Rnw:82-84
###################################################
ex.acgh <-
    aCGH.process(ex.acgh, chrom.remove.threshold = 23, prop.missing = .25, sample.quality.threshold = .4, unmapScreen=TRUE, dupRemove = FALSE)


###################################################
### code chunk number 3: aCGH.Rnw:89-90
###################################################
log2.ratios.imputed(ex.acgh) <- impute.lowess(ex.acgh, maxChrom=24)


###################################################
### code chunk number 4: aCGH.Rnw:96-99
###################################################
data(colorectal)
colorectal
summary(colorectal)


###################################################
### code chunk number 5: aCGH.Rnw:104-105
###################################################
plot(colorectal)


###################################################
### code chunk number 6: aCGH.Rnw:113-115
###################################################
sample.names(colorectal)
phenotype(colorectal)[1:4,]


###################################################
### code chunk number 7: aCGH.Rnw:124-132
###################################################
datadir <- system.file("examples", package = "aCGH")
latest.mapping.file <-
	file.path(datadir, "human.clones.info.Jul03.txt")
ex.acgh <-
	aCGH.read.Sprocs(dir(path = datadir,pattern = "sproc",
			full.names = TRUE), latest.mapping.file,
			chrom.remove.threshold = 23)
ex.acgh


###################################################
### code chunk number 8: aCGH.Rnw:140-141
###################################################
plot(ex.acgh)


###################################################
### code chunk number 9: aCGH.Rnw:151-152
###################################################
cr <- colorectal[ ,1:3]


###################################################
### code chunk number 10: aCGH.Rnw:165-166
###################################################
plotGenome(ex.acgh, samples=2, Y = FALSE)


###################################################
### code chunk number 11: aCGH.Rnw:184-211
###################################################
## Determining hmm states of the clones. In the interest of time, 
##we have commented this step out and used pre-computed results. 

##hmm(ex.acgh) <- find.hmm.states(ex.acgh)
hmm(ex.acgh) <- ex.acgh.hmm

## Merging hmm states

hmm.merged(ex.acgh) <-
   mergeHmmStates(ex.acgh, model.use = 1, minDiff = .25)

## Calculating the standard deviations for each array. Standard error is 
##calculated for each region and then averaged across regions. The final 
##SDs for each samples are contained in sd.samples(exa.acgh)$madGenome.

sd.samples(ex.acgh) <- computeSD.Samples(ex.acgh)

## Finding the genomic events associated with each sample using 
##results of the partitioning into the states.

genomic.events(ex.acgh) <- find.genomic.events(ex.acgh)

## Plotting and printing the hmm states either to the screen or into the 
##postscript file. Each chromosome for each sample is plotted on a separate
##page

##postscript("hmm.states.temp.ps");plotHmmStates(ex.acgh, sample.ind=1);dev.off()


###################################################
### code chunk number 12: aCGH.Rnw:217-218
###################################################
plotHmmStates(colorectal, sample.ind = 1, chr = 1)


###################################################
### code chunk number 13: aCGH.Rnw:247-248
###################################################
plotFreqStat(colorectal, all = T)


###################################################
### code chunk number 14: aCGH.Rnw:261-262
###################################################
summarize.clones(colorectal)[1:10 ,]


###################################################
### code chunk number 15: aCGH.Rnw:268-274
###################################################
factor <- 3
tbl <- threshold.func(log2.ratios(colorectal),
            posThres=factor*(sd.samples(colorectal)$madGenome))
rownames(tbl) <- clone.names(colorectal)
colnames(tbl) <- sample.names(colorectal)
tbl[1:5,1:5]


###################################################
### code chunk number 16: aCGH.Rnw:279-281
###################################################
col.fga <- fga.func(colorectal, factor=3,chrominfo=human.chrom.info.Jul03)
cbind(gainP=col.fga$gainP,lossP=col.fga$lossP)[1:5,]


###################################################
### code chunk number 17: aCGH.Rnw:297-304
###################################################
colnames(phenotype(colorectal))
sex <- phenotype(colorectal)$sex
sex.na <- !is.na(sex)
index.clones.use <- which(clones.info(colorectal)$Chrom < 23)
colorectal.na <- colorectal[ index.clones.use,sex.na , keep=TRUE]
dat <- log2.ratios.imputed(colorectal.na)
resT.sex <- mt.maxT(dat, sex[sex.na], test = "t.equalvar", B = 1000)


###################################################
### code chunk number 18: aCGH.Rnw:310-312
###################################################
plotFreqStat(colorectal.na, resT.sex, sex[sex.na], factor=3, titles =
             c("Female", "Male"), X = FALSE, Y = FALSE)


###################################################
### code chunk number 19: aCGH.Rnw:322-325
###################################################
plotSummaryProfile(colorectal, response = sex,
                   titles = c("Female", "Male"),
                   X = FALSE, Y = FALSE, maxChrom = 22)


###################################################
### code chunk number 20: aCGH.Rnw:336-344
###################################################
factor <- 3
minChanged <- 0.1
gainloss <- gainLoss(log2.ratios(colorectal)[,sex.na], cols=1:length(which(sex.na)), thres=(factor*(sd.samples(colorectal)$madGenome))[sex.na])
ind.clones.use <- which(gainloss$gainP >= minChanged | gainloss$lossP>= minChanged & clones.info(colorectal)$Chrom < 23)
colorectal.na <- colorectal[ind.clones.use,sex.na, keep=TRUE]
dat <- log2.ratios.imputed(colorectal.na)
resT.sex <- mt.maxT(dat, sex[sex.na],test = "t.equalvar", B = 1000)



###################################################
### code chunk number 21: aCGH.Rnw:350-351
###################################################
plotFreqStat(colorectal.na, resT.sex, sex[sex.na],factor=factor,titles = c("Male", "Female"))


###################################################
### code chunk number 22: aCGH.Rnw:361-369
###################################################
time <- rexp(ncol(colorectal), rate = 1 / 12)
events <- rbinom(ncol(colorectal), size = 1, prob = .5)
surv.obj <- Surv(time, 	events)
surv.obj
stat.coxph <-
  aCGH.test(colorectal, surv.obj, test = "coxph",
	p.adjust.method = "fdr")
stat.coxph[1:10 ,]


###################################################
### code chunk number 23: aCGH.Rnw:374-376
###################################################
plotFreqStat(colorectal, stat.coxph, events, titles =
             c("Survived/Censored", "Dead"), X = FALSE, Y = FALSE)


###################################################
### code chunk number 24: aCGH.Rnw:386-394
###################################################
age <- phenotype(colorectal)$age
age.na <- which(!is.na(age))
age <- age[age.na]
colorectal.na <- colorectal[, age.na]
stat.age <-
  aCGH.test(colorectal.na, age, test = "linear.regression",
            p.adjust.method = "fdr")
stat.age[1:10 ,]


###################################################
### code chunk number 25: aCGH.Rnw:399-401
###################################################
plotFreqStat(colorectal.na, stat.age, ifelse(age < 70, 0, 1), titles =
             c("Young", "Old"), X = FALSE, Y = FALSE)


###################################################
### code chunk number 26: aCGH.Rnw:409-418
###################################################
sex <- phenotype(colorectal)$sex
sex.na <- !is.na(sex)
index.clones.use <- which(clones.info(colorectal.na)$Chrom < 23)
colorectal.na <- colorectal[ index.clones.use,sex.na , keep=TRUE]
dat <- log2.ratios.imputed(colorectal.na)
resT.sex <- mt.maxT(dat, sex[sex.na], test = "t.equalvar", B = 1000)

sex.tbl <- summarize.clones(colorectal.na, resT.sex, sex[sex.na], titles = c("Male", "Female"))
sex.tbl[1:5,]


###################################################
### code chunk number 27: aCGH.Rnw:435-445
###################################################
par(mfrow=c(2,1))
clusterGenome(colorectal.na, response = sex[sex.na],
		titles = c("Female", "Male"),
		byclass = FALSE, showaber = TRUE, vecchrom = c(4,8,9),
		dendPlot = FALSE, imp = FALSE)
clusterGenome(colorectal.na, response = sex[sex.na],
		titles = c("Female", "Male"),
		byclass = TRUE, showaber = TRUE, vecchrom = c(4,8,9),
		dendPlot = FALSE, imp = FALSE)



