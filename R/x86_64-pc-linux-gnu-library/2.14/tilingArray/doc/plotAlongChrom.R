### R code from vignette source 'plotAlongChrom.Rnw'

###################################################
### code chunk number 1: load
###################################################
library("grid")
library("RColorBrewer")
library("tilingArray")


###################################################
### code chunk number 2: errorReporting
###################################################
options(error=recover, warn=0, digits=3)


###################################################
### code chunk number 3: showTilingArrayData
###################################################
data("segnf")
class(segnf)
ls(segnf)
segnf$"1.+"
head(segnf$"1.+"@y)
dim(segnf$"1.+"@y)
head(segnf$"1.+"@x)
length(segnf$"1.+"@x)
segnf$"1.+"@logLik
segnf$"1.+"@nrSegments
head(segnf$"1.+"@breakpoints[[segnf$"1.+"@nrSegments]])


###################################################
### code chunk number 4: gffSub
###################################################
data(gffSub)
head(gffSub)


###################################################
### code chunk number 5: alongChromDot1
###################################################
grid.newpage()
plotAlongChrom(segnf,chr=1, coord=c(35000,50000),what="dots", gff=gffSub)


###################################################
### code chunk number 6: alongChromDot2
###################################################
segObj = new.env(parent = baseenv())
nmLabel = colnames(segnf$"1.+"@y)
lab = gsub("\\d","",nmLabel)
for(nm in paste(1,c("+","-"),sep=".")){
    s = get(nm,env = segnf)
    rpY =  tapply(1:length(lab),lab,function(i)rowMeans(s@y[,i]))
    s@y = do.call(cbind,rpY)
    assign(nm,s,segObj)
}
grid.newpage()
plotAlongChrom(segObj,chr=1, coord=c(35000,50000),what="dots", gff=gffSub,sepPlot = T)


###################################################
### code chunk number 7: alongChromHeatmap1
###################################################
grid.newpage() 
plotAlongChrom(segnf,chr=1, coord=c(35000,50000),what="heatmap", gff=gffSub,
         rowNamesHeatmap=nmLabel,makeRasterImage=FALSE)


###################################################
### code chunk number 8: alongChromHeatmap2
###################################################
grid.newpage() 
plotAlongChrom(segnf,chr=1, coord=c(35000,50000),what="heatmap", gff=gffSub,
         rowNamesHeatmap=nmLabel,makeRasterImage=TRUE)


###################################################
### code chunk number 9: alongChromHeatmap3
###################################################
grid.newpage() 
plotAlongChrom(segnf,chr=1, coord=c(35000,50000),what="heatmap", gff=gffSub,
         rowNamesHeatmap=nmLabel,makeRasterImage=TRUE,
         colHeatmap = colorRamp(brewer.pal(9, "Blues")))


