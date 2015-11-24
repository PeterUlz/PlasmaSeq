library(tilingArray)
options(error=recover)

for(f in dir("~/madman/Rpacks/tilingArray/R", pattern=".R$", full.names=TRUE))source(f)



set.seed(123)

## TEST 1: old versus new
if(!TRUE){
  lseg  = 12
  steps = seq(1, 4, length=20)
  steps = steps * sign(runif(length(steps))-0.5)
  dat   = cbind(rep(cumsum(steps), each=lseg) + rnorm(lseg*length(steps)))
  
  par(mfrow=c(2,1))
  
  seg1 = segment(dat, maxseg=30, maxk=nrow(dat))
  seg2 = confint(seg1, 20)
  
  plot(seg2, 20, pch=".")
  
  old1 = findSegments(dat, maxcp=30, maxk=nrow(dat))
  old2 = confint.segmentation(old1, 20)
  plot.segmentation(old2, 20, pch=".")
}

## TEST 2: dependence on number of replicates
lseg  = 20
steps = c(1, 1.05)
steps = steps * sign(runif(length(steps))-0.5)
nrrep = 30
res = matrix(NA, ncol=3, nrow=nrrep)
for(r in 1:nrrep) {
  dat   = sapply(1:r, function(i) rep(cumsum(steps), each=lseg) + rnorm(lseg*length(steps)))
  s = segment(dat, maxseg=5, maxk=nrow(dat))
  s = confint(s, 2)
  res[r, ] = s@breakpoints[[2]]
}
matplot(res, type="b", pch=16)
