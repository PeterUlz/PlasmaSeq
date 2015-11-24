args <- commandArgs()

dataTable <-read.table(args[4], header=TRUE);

ratio<-data.frame(dataTable)

for (i in c(1:24)) {
        png(filename = paste(args[4],".", i, ".segmented.png",sep = ""), width = 1180, height = 1180, units = "px", pointsize = 20, bg = "white", res = NA)
	tt <- which(ratio$Chromosome==i)
	if (length(tt)>0) {
	 plot(ratio$Position[tt],ratio$LogRatio[tt],ylim = c(-2,2),xlab = paste ("position, chr",i),ylab = "log2-ratios",pch = ".",col = colors()[220], cex=3)
         points(ratio$Position[tt],ratio$Summary[tt], pch = ".", col = colors()[88],cex=5)
	 tt <- which(ratio$Chromosome==i  & ratio$Summary > 0.2 )
	 points(ratio$Position[tt],ratio$Summary[tt], pch = ".", col = colors()[136],cex=5)
	 tt <- which(ratio$Chromosome==i  & ratio$Summary < -0.2)
	 points(ratio$Position[tt],ratio$Summary[tt], pch = ".", col = colors()[461],cex=5)
	}
}

dev.off()
warnings()

