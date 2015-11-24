args <- commandArgs()

dataTable <-read.table(args[4], header=TRUE);

ratio<-data.frame(dataTable)

png(filename = paste(args[4],".segments.png",sep = ""), width = 1180, height = 1180,
    units = "px", pointsize = 20, bg = "white", res = NA)
plot(1:10)
op <- par(mfrow = c(5,5))
for (i in c(1:24)) {
	tt <- which(ratio$Chromosome==i)
	if (length(tt)>0) {
	 plot(ratio$Position[tt],ratio$LogRatio[tt],ylim = c(-2,2),xlab = paste ("position, chr",i),ylab = "log2-ratios",pch = ".",col = colors()[200])
         points(ratio$Position[tt],ratio$Summary[tt], pch = ".", col = colors()[88],cex=1)
	 tt <- which(ratio$Chromosome==i  & ratio$Summary > 0.2 )
	 points(ratio$Position[tt],ratio$Summary[tt], pch = ".", col = colors()[136],cex=1)
	 tt <- which(ratio$Chromosome==i  & ratio$Summary < -0.2)
	 points(ratio$Position[tt],ratio$Summary[tt], pch = ".", col = colors()[461],cex=1)
	}
}

dev.off()


