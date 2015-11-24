
## if(dev.cur() <= 1) get(getOption("device"))()

opar <-
    par(ask = interactive() &&
        (.Device %in% c("X11", "GTK", "windows","quartz")))

## Reading Sproc files

datadir <- system.file("examples", package = "aCGH")
latest.mapping.file <-
      file.path(datadir, "human.clones.info.Jul03.txt")
ex.acgh <-
	aCGH.read.Sprocs(dir(path = datadir,pattern = "sproc",
			full.names = TRUE), latest.mapping.file,
			chrom.remove.threshold = 23)
ex.acgh

## Load the colorectal example

data(colorectal)
colorectal

## Basic heatmap plot for batch of aCGH Sproc files

plot(colorectal)

## Subsetting example

colorectal[1:1000 ,1:3]

## Basic plot for the ordered log2 ratios along the genome
## Plot samples in the order of descending quality 

order.quality <- order(sd.samples(colorectal)$madGenome)
par(mfrow = c(2, 1))
pdf("plotGenome.orderByQuality.pdf")
for (i in order.quality)
    plotGenome(colorectal, samples = i, Y = FALSE)
dev.off()

#cluster all samples using imputed data on all chromosomes (autosomes
## and X):

par(mfrow = c(1, 1))
clusterGenome(colorectal, dendPlot = FALSE)

## Plotting the hmm states
plotHmmStates(colorectal, 1)

## Plotting summary of the tumor profiles

sex <- phenotype(colorectal)$sex
plotSummaryProfile(colorectal, response = sex,
                   titles = c("Male", "Female"))

## Testing association of clones with sex, followed by
## plotting the p.values 
## Use mt.maxT function from multtest package to test
## differences in group means for each clone grouped by sex
## Plot the result along the genome

sex <- phenotype(colorectal)$sex
sex.na <- !is.na(sex)
colorectal.na <- colorectal[ ,sex.na, keep = TRUE ]
dat <- log2.ratios.imputed(colorectal.na)
resT.sex <- mt.maxT(dat, sex[sex.na], test = "t", B = 1000)
## WARNING: This takes some time to plot
plotFreqStat(colorectal.na, resT.sex, sex[sex.na],
             titles = c("Male", "Female"))

## Adjust the p.values from previous exercise with "fdr"
## method and plot them
resT.sex.fdr <- resT.sex
resT.sex.fdr$adjp <- p.adjust(resT.sex.fdr$rawp, "fdr")
plotFreqStat(colorectal.na, resT.sex.fdr, sex[sex.na],
             titles = c("Male", "Female"), X = FALSE)

## Derive statistics and p-values for testing the linear association of
## age with the log2 ratios of each clone along the tumors

##age <- phenotype(colorectal)$age
##age.na <- !is.na(age)
##colorectal.na <- colorectal[ ,age.na, keep = T ]
##dat <- log2.ratios.imputed(colorectal.na)
##stat.age <-
##    sapply(1:nrow(dat),
##           function(i) {
               
##               if (i %% 100 == 0)
##                   cat(i, "\n")
##               lm.fit <-
##                   summary(lm(dat[i,] ~ age[age.na]))
##               c(lm.fit$fstatistic[1],
##                 1 - pf(lm.fit$fstatistic[1],
##                        lm.fit$fstatistic[2],
##                        lm.fit$fstatistic[3])
##                 )
               
##           }
##           )

#### Make resT
##resT.age <- data.frame(index = 1:ncol(stat.age),
##                       teststat = stat.age[1,],
##                       rawp = stat.age[2,],
##                       adjp = p.adjust(stat.age[2,], "fdr"))
##resT.age <- resT.age[order(resT.age$adjp),]

#### WARNING: This takes some time to plot
##plotFreqStat(colorectal.na, resT.age, age[age.na])

par(opar)
