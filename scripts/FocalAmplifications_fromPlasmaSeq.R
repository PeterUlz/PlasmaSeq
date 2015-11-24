args <- commandArgs()

#Usage cat FocalAmplifications_fromPlasmaSeq.R | R --slave --args <segmented> <segments> <Type [Prostate | Breast]> <genes> 
ratio<-read.table(args[4], header=TRUE);
segments<-read.table(args[5], header=TRUE);
output_prefix<-args[6]
type<-args[7]
gene_count_cutoff<-strtoi(args[8])
ref_prefix<-args[9]

genes_file <-read.table(paste(ref_prefix,"/genes_unique.txt",sep=""), header=TRUE);

if (type == "Prostate") {
  genes_amp <-read.table(paste(ref_prefix,"/prostate_amplification_genes.txt",sep=""), header=TRUE, stringsAsFactors=FALSE);
  genes_del <-read.table(paste(ref_prefix,"/prostate_deletion_genes.txt",sep=""), header=TRUE, stringsAsFactors=FALSE);
}
if (type == "Breast") {
  genes_amp <-read.table(paste(ref_prefix,"/breast_amplification_genes.txt",sep=""), header=TRUE, stringsAsFactors=FALSE);
  genes_del <-read.table(paste(ref_prefix,"/breast_deletion_genes.txt",sep=""), header=TRUE, stringsAsFactors=FALSE);
}

dgv_variants<-read.table(paste(ref_prefix,"/dgv_variants_hg19_161014.txt",sep=""), header=TRUE, stringsAsFactors=FALSE);
segdup <-read.table(paste(ref_prefix,"/segmental_duplications_merged.bed",sep=""), header=TRUE);

chr_arms <- read.table(paste(ref_prefix,"/chrom_arm_hg19.txt",sep=""), header = TRUE)


cat(paste(type,"cancer","\n"))
cat(paste("Using Gene count cutoff:",toString(gene_count_cutoff),"\n", sep="\t"))

#-------------------------------------------------------------------------------------------------------------------
#mark segments covering genes from important genes in putative amplifications.
#make list entries for each gene
genes_amp_names<-rep("NA", nrow(segments))
genes_amp_list<-rep("NA", nrow(segments))
genes_namelist<-rep("NA", nrow(segments))
segments<-cbind(segments, genes_amp_names, stringsAsFactors=FALSE)
segments<-cbind(segments, genes_amp_list, stringsAsFactors=FALSE)
segments<-cbind(segments, genes_namelist, stringsAsFactors=FALSE)
for (i in 1:nrow(segments)){
   # full length in segment only
   count<-which(genes_amp$chrom == segments$Chromosome[i] & genes_amp$cdsStart > segments$Start[i] &  genes_amp$cdsEnd < segments$End[i] )
   if (length(count) > 0) {
      #levels(segments$genes_amp_names)<-c(levels(segments$genes_amp_names),paste(genes_amp$name2[count], sep="", collapse=','))
      segments$genes_amp_names[i]<-paste(genes_amp$name2[count], sep="", collapse=',')
      #levels(segments$genes_amp_list)<-c(levels(segments$genes_amp_list),paste(genes_amp$list[count], sep="", collapse=','))
      segments$genes_amp_list[i]<-paste(genes_amp$list[count], sep="", collapse=',')
      segments$genes_namelist[i]<-paste(genes_amp$name2[count],"(",genes_amp$list[count],")", sep="", collapse=',')
   }
   # annotate genes with partial overlap
   count<-which(genes_amp$chrom == segments$Chromosome[i] & ((genes_amp$cdsStart < segments$Start[i] &  genes_amp$cdsEnd > segments$Start[i]) | (genes_amp$cdsStart < segments$End[i] &  genes_amp$cdsEnd > segments$End[i])) )
   if (length(count) > 0) {
      #levels(segments$genes_amp_names)<-c(levels(segments$genes_amp_names),paste(genes_amp$name2[count], sep="", collapse=','))
      partial_gene_names<-paste(genes_amp$name2[count], "(partial)", sep="")
      segments$genes_amp_names[i]<-paste(c(segments$genes_amp_names[i],partial_gene_names), sep="", collapse=',')
      #levels(segments$genes_amp_list)<-c(levels(segments$genes_amp_list),paste(genes_amp$list[count], sep="", collapse=','))
      segments$genes_amp_list[i]<-paste(genes_amp$list[count], sep="", collapse=',')
      segments$genes_namelist[i]<-paste(genes_amp$name2[count],"(",genes_amp$list[count],")(partial)", sep="", collapse=',')
   }
}
#-------------------------------------------------------------------------------------------------------------------
#mark segments covering genes from important genes in putative deletions.
genes_del_names<-rep("NA", nrow(segments))
segments<-cbind(segments, genes_del_names, stringsAsFactors=FALSE)
for (i in 1:nrow(segments)){
  count<-which(genes_del$chrom == segments$Chromosome[i] & genes_del$cdsStart > segments$Start[i] &  genes_del$cdsEnd < segments$End[i] )
  if (length(count) > 0) {
     #levels(segments$genes_del_names)<-c(levels(segments$genes_del_names),paste(genes_del$name2[count], sep="", collapse=','))
     segments$genes_del_names[i]<-paste(genes_del$name2[count], sep="", collapse=',')
  }
   # annotate genes with partial overlap
   count<-which(genes_del$chrom == segments$Chromosome[i] & ((genes_del$cdsStart < segments$Start[i] &  genes_del$cdsEnd > segments$Start[i]) | (genes_del$cdsStart < segments$End[i] &  genes_del$cdsEnd > segments$End[i])) )
   if (length(count) > 0) {
      #levels(segments$genes_del_names)<-c(levels(segments$genes_del_names),paste(genes_del$name2[count], sep="", collapse=','))
      partial_gene_names<-paste(genes_del$name2[count], "(partial)", sep="")
      segments$genes_del_names[i]<-paste(c(segments$genes_del_names[i],partial_gene_names), sep="", collapse=',')
   }
}

#-------------------------------------------------------------------------------------------------------------------
#mark segments covering genes (only for segments smaller than 20Mbp)
genes<-rep("NA", nrow(segments))
gene_count<-rep(0, nrow(segments))
segments<-cbind(segments, genes, gene_count)
for (i in 1:nrow(segments)){
   #only compute gene content for amplified segments smaller than 20Mb
   if ( (segments$Log2Ratio[i] > 0.2 | segments$Log2Ratio[i] < -0.2) &  (segments$End[i]-segments$Start[i])<20000000)
   {
      count<-which(genes_file$chrom == segments$Chromosome[i] & genes_file$txStart > segments$Start[i] &  genes_file$txEnd < segments$End[i] )
      if (length(count) > 0) {
         levels(segments$genes)<-c(levels(segments$genes),paste(genes_file$name2[count], sep="", collapse=','))
         segments$genes[i]<-paste(genes_file$name2[count], sep="", collapse=',')
	 segments$gene_count[i]<-length(count)
      }
   }
}

#-------------------------------------------------------------------------------------------------------------------
#mark segments covering segmental duplications with >50%
percent_segdup<-rep(0, nrow(segments))
segments<-cbind(segments, percent_segdup)
for (i in 1:nrow(segments)){
   #only compute segdup content for amplified segments smaller than 20Mb
   cum_length<-0
   if (  (segments$Log2Ratio[i] > 0.2 | segments$Log2Ratio[i] < -0.2) & (segments$End[i]-segments$Start[i])<20000000)
   {
      tt<-which(segdup$Chrom == segments$Chromosome[i] & segdup$Start > segments$Start[i] &  segdup$End < segments$End[i] )
      if (length(tt) > 0) {
         cum_length<-sum(segdup$End[tt] - segdup$Start[tt])   
      }
      tt<-which(segdup$Chrom == segments$Chromosome[i] & segdup$Start < segments$Start[i] &  segdup$End < segments$End[i] & segdup$End > segments$Start[i])
      if (length(tt) > 0) {
	 cum_length<-cum_length+sum(segdup$End[tt] - segments$Start[i] )   
      }
      tt<-which(segdup$Chrom == segments$Chromosome[i] & segdup$Start > segments$Start[i] &  segdup$End > segments$End[i] & segdup$Start < segments$End[i])
      if (length(tt) > 0) {
         cum_length<-cum_length+sum( segments$End[i] - segdup$Start[tt])
      }
      tt<-which(segdup$Chrom == segments$Chromosome[i] & segdup$Start < segments$Start[i] &  segdup$End > segments$End[i])
      if (length(tt) > 0) {
         cum_length<-cum_length+sum( segments$End[i] - segments$Start[i] )
      }
   }
   cum_length
   segments$percent_segdup[i]<-cum_length / (segments$End[i]-segments$Start[i])
}

#-------------------------------------------------------------------------------------------------------------------
#mark segments covering variants from dgv when variant is within 500kbp of start and end of segment
dgv<-rep("NA", nrow(segments))
segments<-cbind(segments, dgv, stringsAsFactors=FALSE)
dgv_type<-rep("NA", nrow(segments))
segments<-cbind(segments, dgv_type, stringsAsFactors=FALSE)

for (i in 1:nrow(segments)){
   #For very small segments be more strict for annotating DGV entries
   if (  (segments$Log2Ratio[i] > 0.2 | segments$Log2Ratio[i] < -0.2) & (segments$End[i]-segments$Start[i])<500000)
   {
      count<-which(dgv_variants$chr == segments$Chromosome[i] & dgv_variants$start > (segments$Start[i]-10000) &  dgv_variants$start < (segments$Start[i]+10000) & dgv_variants$end > (segments$End[i]-10000) &  dgv_variants$end < (segments$End[i]+10000) )
      #if more than 10 entries do not print individual accession numbers
      if (length(count) > 10){
         segments$dgv[i]<-">10"
         segments$dgv_type[i]<-">10"
      }
      else if (length(count) > 0) {
         segments$dgv[i]<-paste(dgv_variants$variantaccession[count],sep="",collapse=",")
         segments$dgv_type[i]<-paste(dgv_variants$variantsubtype[count],sep="",collapse=",")
      }
   }
   #only compute cnv content for amplified segments smaller than 20Mb
   else if (  (segments$Log2Ratio[i] > 0.2 | segments$Log2Ratio[i] < -0.2) & (segments$End[i]-segments$Start[i])<20000000)
   {
      count<-which(dgv_variants$chr == segments$Chromosome[i] & dgv_variants$start > (segments$Start[i]-100000) &  dgv_variants$start < (segments$Start[i]+100000) & dgv_variants$end > (segments$End[i]-100000) &  dgv_variants$end < (segments$End[i]+100000) )
      #if more than 10 entries do not print individual accession numbers
      if (length(count) > 10){
         segments$dgv[i]<-">10"
         segments$dgv_type[i]<-">10"
      }
      else if (length(count) > 0) {
         segments$dgv[i]<-paste(dgv_variants$variantaccession[count],sep="",collapse=",")
         segments$dgv_type[i]<-paste(dgv_variants$variantsubtype[count],sep="",collapse=",")
      }
   }
}

#-------------------------------------------------------------------------------------------------------------------
#compute chr arm average and add columns "Chromosome arm average Log2" and "Chr. Arm. adjusted Log2Ratio" to data frame segments
log2_average<-rep(0, nrow(chr_arms))
chr_arm_average<-rep(0, nrow(segments))
Log2Ratio_chr_arm_adjust<-rep(0, nrow(segments))
chr_arms<-cbind(chr_arms, log2_average)
segments<-cbind(segments, chr_arm_average)
segments<-cbind(segments, Log2Ratio_chr_arm_adjust)
for (i in 1:nrow(chr_arms)){
   tt<-which(segments$Chromosome == chr_arms$Chrom[i] & segments$Start >= chr_arms$Start[i] & segments$End <= chr_arms$End[i] )
   if (length(tt) > 0){
      log2_ratio_sum<-sum((segments$End[tt]-segments$Start[tt])*segments$Log2Ratio[tt])
      length_arm<-max(segments$End[tt])-min(segments$Start[tt])
      chr_arms$log2_average[i]<-log2_ratio_sum / length_arm
      segments$chr_arm_average[tt]<-log2_ratio_sum / length_arm
      segments$Log2Ratio_chr_arm_adjust[tt]<-segments$Log2Ratio[tt] - segments$chr_arm_average[tt] 
   }
}
#-------------------------------------------------------------------------------------------------------------------
#compute log2 average of 20Mbp upstream
log2_average_upstream<-rep(0, nrow(segments))
Log2Ratio_20Mbp_upstream_adjust<-rep(0, nrow(segments))
segments<-cbind(segments, log2_average_upstream)
segments<-cbind(segments, Log2Ratio_20Mbp_upstream_adjust)
for (i in 1:nrow(segments)){
   log2_ratio_sum<-0
   #take every segment completely inside 20Mbp
   tt_inside<-which(segments$Chromosome == segments$Chromosome[i] & segments$Start < segments$Start[i] & segments$End < segments$End[i] & segments$Start > (segments$Start[i]-20000000))
   if (length(tt_inside) > 0){
      log2_ratio_sum<-sum((segments$End[tt_inside]-segments$Start[tt_inside])*segments$Log2Ratio[tt_inside])
   }
   #take partial segments between 20Mbp and beyond
   tt_border<-which(segments$Chromosome == segments$Chromosome[i] & segments$Start < segments$Start[i] & segments$Start < (segments$Start[i]-20000000) & segments$End > (segments$Start[i]-20000000))
   if (length(tt_border) > 0){
      log2_ratio_sum<-log2_ratio_sum + ( segments$End[tt_border] - (segments$Start[i]-20000000))*segments$Log2Ratio[tt_border]
      segments$log2_average_upstream[i]<-log2_ratio_sum / 20000000
      segments$Log2Ratio_20Mbp_upstream_adjust[i]<-segments$Log2Ratio[i] - segments$log2_average_upstream[i]
   }
   #if segment in question is closer than 20Mb at the chromosome edge just count internal segments
   if (length(tt_inside) != 0 & length(tt_border) == 0){
      segments$Log2Ratio_20Mbp_upstream_adjust[i]<-segments$Log2Ratio[i]
      log2_ratio_sum<-sum((segments$End[tt_inside]-segments$Start[tt_inside])*segments$Log2Ratio[tt_inside])
      length_segment<-max(segments$End[tt_inside]) - min(segments$Start[tt_inside])
      segments$log2_average_upstream[i]<-log2_ratio_sum / length_segment
      segments$Log2Ratio_20Mbp_upstream_adjust[i]<-segments$Log2Ratio[i] - segments$log2_average_upstream[i]
   }
   if (length(tt_inside) == 0 & length(tt_border) == 0){
      segments$Log2Ratio_20Mbp_upstream_adjust[i]<-segments$Log2Ratio[i]
   }
}
#-------------------------------------------------------------------------------------------------------------------
#compute log2 average of 20Mbp downstream
log2_average_downstream<-rep(0, nrow(segments))
Log2Ratio_20Mbp_downstream_adjust<-rep(0, nrow(segments))
segments<-cbind(segments, log2_average_downstream)
segments<-cbind(segments, Log2Ratio_20Mbp_downstream_adjust)
for (i in 1:nrow(segments)){
   log2_ratio_sum<-0
   #take every segment completely inside 20Mbp
   tt_inside<-which(segments$Chromosome == segments$Chromosome[i] & segments$Start > segments$Start[i] & segments$End > segments$End[i] & segments$End < (segments$End[i]+20000000))
   if (length(tt_inside) > 0){
      log2_ratio_sum<-sum((segments$End[tt_inside]-segments$Start[tt_inside])*segments$Log2Ratio[tt_inside])
   }
   #take partial segments between 20Mbp and beyond
   tt_border<-which(segments$Chromosome == segments$Chromosome[i] & segments$Start > segments$Start[i] & segments$End > (segments$End[i]+20000000) & segments$Start < (segments$End[i]+20000000))
   if (length(tt_border) > 0){
      log2_ratio_sum<-log2_ratio_sum + ((segments$End[i]+20000000) - segments$Start[tt_border])*segments$Log2Ratio[tt_border]
      segments$log2_average_downstream[i]<-log2_ratio_sum / 20000000
      segments$Log2Ratio_20Mbp_downstream_adjust[i]<-segments$Log2Ratio[i] - segments$log2_average_downstream[i]
   }
   #if segment in question is closer than 20Mb at the chromosome edge just count internal segments
   if (length(tt_inside) != 0 & length(tt_border) == 0){
      segments$Log2Ratio_20Mbp_downstream_adjust[i]<-segments$Log2Ratio[i]
      log2_ratio_sum<-sum((segments$End[tt_inside]-segments$Start[tt_inside])*segments$Log2Ratio[tt_inside])
      length_segment<-max(segments$End[tt_inside]) - min(segments$Start[tt_inside])
      segments$log2_average_downstream[i]<-log2_ratio_sum / length_segment
      segments$Log2Ratio_20Mbp_downstream_adjust[i]<-segments$Log2Ratio[i] - segments$log2_average_downstream[i]
   }
   if (length(tt_inside) == 0 & length(tt_border) == 0){
      segments$Log2Ratio_20Mbp_downstream_adjust[i]<-segments$Log2Ratio[i]
   }
}

#-------------------------------------------------------------------------------------------------------------------
#draw pictures and write data to output file
write.table(t(c("Chromosome", "Start","End","Log2Ratio","Genelist", "Gene","Gene (Genelist)","Chrom Arm average", "Log2Ratio chrom arm adjusted", "20Mbp upstream average", "Log2 20Mb upstream adjusted", "20Mbp downstream average", "Log2 20Mb downstream adjusted", "Other genes", "Gene count", "Percent SegDup", "DGV accessions", "DGV types")), file=paste(output_prefix, "_focal_only.csv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t" )
for (i in c(1:24)) {
        png(filename = paste(output_prefix,".", i, ".segmented.png",sep = ""), width = 1180, height = 1180, units = "px", pointsize = 20, bg = "white", res = NA)
	tt <- which(ratio$Chromosome==i)
	if (i == 23){
           i = "X"
        }
	if (i == 24){
           i = "Y"
        }

	if (length(tt)>0) {
	 plot(ratio$Position[tt],ratio$LogRatio[tt],ylim = c(-2,2),xlab = paste ("position, chr",i),ylab = "log2-ratios",pch = ".",col = colors()[220], cex=3)
	 tt<-which(segments$Chromosome == paste("chr", i, sep=""))
	 segments(segments$Start[tt]+1, segments$Log2Ratio[tt], segments$End[tt], segments$Log2Ratio[tt], col=colors()[88], lwd=5)
	 tt<-which(segments$Chromosome == paste("chr", i, sep="") & segments$Log2Ratio > 0.2)
	 segments(segments$Start[tt]+1, segments$Log2Ratio[tt], segments$End[tt], segments$Log2Ratio[tt], col=colors()[136], lwd=5)
         tt<-which(segments$Chromosome == paste("chr", i, sep="") & segments$Log2Ratio < -0.2)
	 segments(segments$Start[tt]+1, segments$Log2Ratio[tt], segments$End[tt], segments$Log2Ratio[tt], col=colors()[461], lwd=5)

	 #only use segments having log-ratio > 0.2; adjusted 20Mb log2ratio >0.2 (upstream and downstream); segments < 20Mb length; segments containing less than (cutoff) genes; less than 50% segmental duplications and no cnvs
         known_count<-which(segments$Chromosome == paste("chr", i, sep="") & segments$Log2Ratio > 0.2  & segments$genes_amp_names != "NA" & (segments$End-segments$Start)<20000000 & segments$Log2Ratio_20Mbp_upstream_adjust > 0.2 & segments$Log2Ratio_20Mbp_downstream_adjust > 0.2  & segments$gene_count < gene_count_cutoff & segments$dgv != ">10" & !grepl("gain",segments$dgv_type) & !grepl("insertion",segments$dgv_type) & !grepl("duplication",segments$dgv_type))
         if (length(known_count)>0) {
            rect(segments$Start[known_count]+1, rep(-1.9, length(known_count)), segments$End[known_count], rep(1.9, length(known_count)), border=NULL, col=rgb(0.9, 0.1, 0.1, 0.3))
            text( (segments$Start[known_count]+segments$End[known_count])/2, rep(2, length(known_count)), labels=segments$genes_amp_names[known_count], cex=0.8)
	         for (y in known_count) {
               write.table(t(c(toString(segments$Chromosome[y]), segments$Start[y], segments$End[y], segments$Log2Ratio[y], segments$genes_amp_list[y], segments$genes_amp_names[y], segments$genes_namelist[y], segments$chr_arm_average[y], segments$Log2Ratio_chr_arm_adjust[y], segments$log2_average_upstream[y], segments$Log2Ratio_20Mbp_upstream_adjust[y], segments$log2_average_upstream[y], segments$Log2Ratio_20Mbp_downstream_adjust[y], toString(segments$genes[y]), segments$gene_count[y], toString(segments$percent_segdup[y]), segments$dgv[y], segments$dgv_type[y])), file=paste(output_prefix, "_focal_only.csv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE,append=TRUE, sep="\t")
            }
         }

	 #only use segments having log-ratio > 0.2; adjusted 20Mb log2ratio >0.58 (upstream and downstream); segments < 20Mb length; segments containing less than (cutoff) genes or less and not already in known genes
         unknown_criteria<-which(segments$Chromosome == paste("chr", i, sep="") & segments$Log2Ratio > 0.2  & segments$genes != "NA" & (segments$End-segments$Start)<20000000 & segments$Log2Ratio_20Mbp_upstream_adjust > 0.58 & segments$Log2Ratio_20Mbp_downstream_adjust > 0.58 & segments$gene_count < gene_count_cutoff & segments$percent_segdup < 0.5 & segments$dgv != ">10" & !grepl("gain",segments$dgv_type) & !grepl("insertion",segments$dgv_type) & !grepl("duplication",segments$dgv_type))
         if (length(known_count) > 0 & length(unknown_criteria) > 0) {
            unknown_count<-setdiff(unknown_criteria,known_count)
         } else {
            unknown_count<-unknown_criteria
         }
         if (length(unknown_count)>0) {
            rect(segments$Start[unknown_count]+1, rep(-1.9, length(unknown_count)), segments$End[unknown_count], rep(1.9, length(unknown_count)), border=NULL, col=rgb(0.9, 0.1, 0.1, 0.3))
            text( (segments$Start[unknown_count]+segments$End[unknown_count])/2, rep(2, length(unknown_count)), labels=segments$genes[unknown_count], cex=0.8)
	         for (y in unknown_count) {
               write.table(t(c(toString(segments$Chromosome[y]), segments$Start[y], segments$End[y], segments$Log2Ratio[y], "Other", toString(segments$genes[y]), "NA",segments$chr_arm_average[y], segments$Log2Ratio_chr_arm_adjust[y], segments$log2_average_upstream[y], segments$Log2Ratio_20Mbp_upstream_adjust[y], segments$log2_average_upstream[y], segments$Log2Ratio_20Mbp_downstream_adjust[y], toString(segments$genes[y]), segments$gene_count[y], toString(segments$percent_segdup[y]), segments$dgv[y], segments$dgv_type[y])), file=paste(output_prefix, "_focal_only.csv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE,append=TRUE, sep="\t")
            }
         }

	 #only use segments having log-ratio < -0.2; adjusted 20Mb log2ratio <-0.2 (upstream and downstream); segments < 20Mb length; segments containing less than (cutoff) genes
         known_del_count<-which(segments$Chromosome == paste("chr", i, sep="") & segments$Log2Ratio < -0.2 & segments$genes_del_names != "NA" & (segments$End-segments$Start)<20000000 & segments$Log2Ratio_20Mbp_upstream_adjust < -0.2 & segments$Log2Ratio_20Mbp_downstream_adjust < -0.2 & segments$gene_count < gene_count_cutoff & segments$percent_segdup < 0.5 & segments$dgv != ">10" & !grepl("loss",segments$dgv_type) & !grepl("deletion",segments$dgv_type))     
         if (length(known_del_count)>0) {
           rect(segments$Start[known_del_count]+1, rep(-1.9, length(known_del_count)), segments$End[known_del_count], rep(1.9, length(known_del_count)), border=NULL, col=rgb(0.1, 0.8, 0.8, 0.3))
           text( (segments$Start[known_del_count]+segments$End[known_del_count])/2, rep(-2, length(known_del_count)), labels=segments$genes_del_names[known_del_count], cex=0.8)
	        for (y in known_del_count) {
               write.table(t(c(toString(segments$Chromosome[y]), segments$Start[y], segments$End[y], segments$Log2Ratio[y], "Deleted Known", toString(segments$genes_del_names[y]),"NA", segments$chr_arm_average[y], segments$Log2Ratio_chr_arm_adjust[y], segments$log2_average_upstream[y], segments$Log2Ratio_20Mbp_upstream_adjust[y], segments$log2_average_upstream[y], segments$Log2Ratio_20Mbp_downstream_adjust[y], toString(segments$genes[y]), segments$gene_count[y], toString(segments$percent_segdup[y]), segments$dgv[y], segments$dgv_type[y])), file=paste(output_prefix, "_focal_only.csv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE,append=TRUE, sep="\t")
           }
         }
   }

}

#Write out all segments plus information
write.table(t(c("Chromosome", "Start","End","Log2Ratio", "Genelist","Gene (Amplified)","Gene (Deleted)","Other Genes","Chrom Arm average", "Log2Ratio chrom arm adjusted", "20Mbp upstream average", "Log2 20Mb upstream adjusted", "20Mbp downstream average", "Log2 20Mb downstream adjusted", "Percent SegDup", "Gene Count", "DGV accessions", "DGV types")), file=paste(output_prefix, "_all.csv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t" )
for (y in 1:nrow(segments)) {
   write.table(t(c(toString(segments$Chromosome[y]), segments$Start[y], segments$End[y], segments$Log2Ratio[y],toString(segments$genes_amp_list[y]), toString(segments$genes_amp_names[y]), toString(segments$genes_del_names[y]),toString(segments$genes[y]),segments$chr_arm_average[y], segments$Log2Ratio_chr_arm_adjust[y], segments$log2_average_upstream[y], segments$Log2Ratio_20Mbp_upstream_adjust[y],segments$log2_average_downstream[y], segments$Log2Ratio_20Mbp_downstream_adjust[y], toString(segments$percent_segdup[y]),toString(segments$gene_count[y]), segments$dgv[y], segments$dgv_type[y])), file=paste(output_prefix, "_all.csv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE,append=TRUE, sep="\t")
}

