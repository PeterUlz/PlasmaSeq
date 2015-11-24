args <- commandArgs()

ratio<-read.table(args[4], header=TRUE);

png(filename = paste(args[4],"_chrscale_blue.png",sep = ""), width = 2280, height = 218,
    units = "px", pointsize = 20, bg = "white", res = NA)
par(mar=c(4,0,0,0))
count<-1
widths<-vector(length=24)
for (i in c(1:22,"X","Y")) {
	tt <- which(ratio$Chromosome==paste("chr",i, sep=""))
	widths[count]<-max(ratio$End[tt])
        count<-count+1
}
widths
nf <- layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24), 1, 24, byrow=TRUE), widths=widths)
for (i in c(1:22,"X","Y")) {
	tt <- which(ratio$Chromosome==paste("chr",i, sep=""))
	if (length(tt)>0) {
	 plot(ratio$Start[tt],ratio$Log2Ratio[tt],ylim = c(-2,2),xlab = paste("chr",i, sep=""),
          ylab = "log2-ratios",pch = ".",col = colors()[88], , yaxt="n",xaxt="n",xlim=c(0,max(ratio$End[tt])))
	 segments(x0=ratio$Start[tt], y0=ratio$Log2Ratio[tt], x1=ratio$End[tt], y1=ratio$Log2Ratio[tt],lwd=5,col = colors()[461], cex=2)
        }
}

dev.off()

