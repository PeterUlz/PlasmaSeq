### R code from vignette source 'seqLogo.Rnw'

###################################################
### code chunk number 1: load
###################################################
library(seqLogo) 


###################################################
### code chunk number 2: makePWM
###################################################
mFile <- system.file("Exfiles/pwm1", package="seqLogo")
m <- read.table(mFile)
m
p <- makePWM(m)


###################################################
### code chunk number 3: slots
###################################################
slotNames(p)
p@pwm
p@ic
p@consensus


###################################################
### code chunk number 4: cosmoArgs
###################################################
args(seqLogo)


###################################################
### code chunk number 5: seqLogo1
###################################################
seqLogo(p)


###################################################
### code chunk number 6: seqLogo2
###################################################
seqLogo(p, ic.scale=FALSE)


