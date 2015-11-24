### R code from vignette source 'findsegments.Rnw'

###################################################
### code chunk number 1: defGenData
###################################################
genData = function(lenx, nrSeg, nrRep=1, stddev=0.1) {
  x  = matrix(as.numeric(NA), nrow=lenx, ncol=nrRep)
  cp = sort(sample(1:floor(lenx/15), nrSeg-1) * 15)
  cpb = c(1, cp, lenx+1)
  s  = 0
  for (j in 2:length(cpb)) {
    sel = cpb[j-1]:(cpb[j]-1)
    s = (.5+runif(1))*sign(rnorm(1))+s
    x[sel, ] = rnorm(length(sel)*nrRep, mean=s, sd=stddev)
  }
  return(list(x=x, cp=cp))
}


###################################################
### code chunk number 2: plotData
###################################################
set.seed(4711)
lenx = 1000
nrSeg = 10
gd = genData(lenx, nrSeg)
plot(gd$x, pch=".")
abline(v=gd$cp, col="blue")


###################################################
### code chunk number 3: loadlib
###################################################
library("tilingArray")


###################################################
### code chunk number 4: findSegments
###################################################
maxseg = 12
maxk = 500
seg  = segment(gd$x, maxk=maxk, maxseg=maxseg)
seg
seg@breakpoints[nrSeg+(-1:1)]
gd$cp


###################################################
### code chunk number 5: checkClaim
###################################################
stopifnot(all(gd$cp==seg@breakpoints[[nrSeg]]))


###################################################
### code chunk number 6: testCP
###################################################
for(i in 1:20){
  gd  = genData(lenx, nrSeg)
  seg = segment(gd$x, maxk=maxk, maxseg=maxseg)
  stopifnot(seg@breakpoints[[nrSeg]][, "estimate"] == gd$cp)
}


###################################################
### code chunk number 7: modelSelect1
###################################################
nrSeg = 22
gd = genData(lenx, nrSeg, nrRep=2, stddev=1/3)
s = segment(gd$x, maxk=lenx, maxseg=as.integer(nrSeg * 2.5))


###################################################
### code chunk number 8: plotSimSeg
###################################################
par(mai=c(1,1,0.1,0.01))
plot(row(gd$x), gd$x, pch=".")


###################################################
### code chunk number 9: penLLsim
###################################################
par(mai=c(1,1,0.1,0.01))
plotPenLL(s)


###################################################
### code chunk number 10: printMaxima
###################################################
which.max(logLik(s, penalty="AIC"))
which.max(logLik(s, penalty="BIC"))


