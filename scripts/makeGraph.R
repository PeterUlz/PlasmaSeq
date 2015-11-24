args <- commandArgs()

dataTable <-read.table(args[4], header=TRUE);

ratio<-data.frame(dataTable)

png(filename = paste(args[4],".png",sep = ""), width = 1180, height = 1180,
    units = "px", pointsize = 20, bg = "white", res = NA)
plot(1:10)
op <- par(mfrow = c(5,5))
for (i in c(1:24)) {
	tt <- which(ratio$Chromosome==i)
	if (length(tt)>0) {
	 plot(ratio$Position[tt],ratio$LogRatio[tt],ylim = c(-2,2),xlab = paste ("position, chr",i),ylab = "log2-ratios",pch = ".",col = colors()[88])
	 tt <- which(ratio$Chromosome==i  & ratio$Summary > 0.2 )
	 points(ratio$Position[tt],ratio$LogRatio[tt],pch = ".",col = colors()[136])
	 tt <- which(ratio$Chromosome==i  & ratio$Summary < -0.2)
	 points(ratio$Position[tt],ratio$LogRatio[tt],pch = ".",col = colors()[461])
	 #tt <- which(ratio$Chromosome==i)
	 #points(ratio$Start[tt],ratio$MedianRatio[tt], pch = ".", col = colors()[24],cex=1)
	}
}

dev.off()

