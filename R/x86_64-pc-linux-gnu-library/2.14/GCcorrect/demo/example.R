
###
#
# Installing
#
###
install.packages(repos=NULL,'~/GCcorrect_0.96.tar.gz')

library("GCcorrect")

#####
# Begin by preparing a chromosome structure.
#
# reference: The reference chromosome,
#            Formats accepted:
#            1. A numeric vector with 4 values:
#               0,1,2,3,NA numeric (representing "A","T","C","G",NA)
#            2. A fasta file-name, or 0,1,2,3,NA numeric (representing "A","T","C","G",NA)
#
# repeats:   Unmappable positions in the chromosome
#            Formats accepted:
#            1. Logical vector in R format marking unmappable positions (0 - mappable, 1- repeat)
#            2. A fasta file-name similar to the reference, with "R" marking unmappable positions
#               (see example in GCcorrect/data/repeat75chr1.fa) 
#            We attach to files in GCcorrect/python/ to help create these files. 
#
#######

# Get the following tar (in unix `curl`):
system('curl http://www.stat.berkeley.edu/share/yuvalb/HCCchr1.tar > HCCchr1.tar')

# Untar and Gunzip to get the following files: 
# chr1.fa.gz  HCC1569_chr1_forward.cln.gz  HCC1569_chr1_reverse.cln.gz  repeat45chr1.fa.gz
system('tar -xf HCCchr1.tar; gunzip HCC/*')
setwd('HCC')
# The README file will give some background on this data...

chr1 = prepareChrom("chr1.RData",reference="chr1.fa",
  repeats = "repeat45chr1.fa",details = list("h_sapiens_asm",1,45))

# Reading the repeats file might give out a warning (file ended unexpectedly)
# That is fine. 

# Each chromosome should have a separate structure. 

#####
# Input reads
#
# A table for the forward strand and a table for reverse strand reads that were mapped to a chromosome. 
# Tables should have two columns:
# 1. 5' end of the read (as mapped to the forward strand)
# 2. Length of the fragment (if available, or 0 otherwise)
#    See comment below on length inputs; 
#    
# We usually keep chromosome as a first column, and cut it away. 
# See GCcorrect/data/chr1.forward and GCcorrect/data/chr1.reverse for the format 
#######

dat1_for = as.matrix(read.table('HCC1569_chr1_forward.cln')[,2:3])
dat1_rev = as.matrix(read.table('HCC1569_chr1_reverse.cln')[,2:3])

####### Comments on Read Input    ##########
# Indexing starts at 1 (not 0).
# We assume index points to the 5' end of the read.
#
# Forward strand this will be the first base of the sequenced read.
# An easy way to check:
# chr1$chrline[dat1_for[1,1] + 0:9] should give first bases of the read.
# 
# Reverse strand index should point to the 5' of the complementary fragment.
# chr1$chrline[dat1_rev[1,1] + 0:9] should give first bases of the reverse complement read. 
# This is not(!) the end of the fragment.
# Some mapping software output this by default. Otherwise just subtract (read_len-1). 
#
# Note: Make sure not to use the same pair-end in both the forward and in the reverse strand.
#       An easy way to do this is to first run the matching of 5' and 3',
#       and then only use the 5' end results and the length information. 
#
# Length  
# We mean the number of bases between 5' end of fragment and the 3' end.
# This is not always the default of the mapping software
# (some exclude the second "read" for technical reasons)
#
# If length information is available, using only forward strand mappings
# would probably be sufficient. For single-stranded data though,
# splitting by strand is important.
#
################

#####
#
# 
# If your machine does not have a lot of memory (>10G),
# you can first try the code on a smaller segment of Chromosome 1:
chr1$isrep = chr1$isrep[1:20000000]
chr1$chrline = chr1$isrep[1:20000000]
chr1$isgc = chr1$isgc[1:20000000]
chr1$map10K = chr1$map10K[1:2000]
chr1$gc10K = chr1$gc10K[1:2000]

dat1_for = dat1_for[dat1_for[,1]<(20000000-1000),]
dat1_rev = dat1_rev[dat1_rev[,1]<(20000000-1000),]
#
#
#####

#
####


# If you are using only forward strand, prepare a dummy variable for the reverse:
# dat1_rev = matrix(1,nc=2,nr=2);
# And make sure to always use strand="F" where relevant

dat1 = prepareReads (dat1_for, dat1_rev,chr_len = length(chr1$chrline), 1,45,'h_sapiens_asm')


# Basic examination of the reads:
basicPlots(dat1,types = c(1,2,3,4))
# Plots types are:
# 1) Density of read counts
# 2) Length density curves - total, forward, reverse
# 3) Counts by genome location
# 4) Comparison of counts on between strands

# Examination of relation between reads and chromosome (mappability, GC)

compDataRefPlots(chr1 ,dat1,types = c(1,2,3,4))
# Plot types are:
# 1) Binned readcounts, mappability and GC across the chromosome
# 2) Plot of counts by mappability
# 4) Plot of GC by mappability
#
# 3) Calculates compares how likely a read is to be mapped in an unmappable location,
#    compared to the proportion of unmappable locations. This is another indication
#    of the quality of this mappability score

# To get reproduceable sample...
set.seed(1000)

# We want to create a manual mask for the genome,
# removing zero stretches and pileups
useSamp = !logical(length(chr1$isgc)) 
zeros_10K = findLongZeros(dat1$for10K,2)
for (i in 1:ncol(zeros_10K)) {
    useSamp[((10000*(zeros_10K[1,i]-1)) + 1) : (10000*(zeros_10K[2,i]))] = F
}
useSamp = useSamp[1:length(chr1$isgc)]

high_10K = which(dat1$for10K>600)
for (i in high_10K) {
  useSamp[10000*(i-1) + 1:10000] = F
}
useSamp = useSamp[1:length(chr1$isgc)]
#####
# Create a sample of chromosome 1
#####
sampCh1 = sampleChrom(chr1,dat1,n=10000000,margin = 1000,len_range = c(), memoryopt = TRUE,useSamp = useSamp) 

########
# Notes:
# Sample sturctures are associated with the chromosomes.
# Conditional mean structures are usually estimated for each chromosome separately,
# but can later be combined together using sumCondMeans
########

######
# We want to find the best GC window for this data
# We stratify locations based on the %GC of the window starting margin bp after the location.
# We compare different window lengths, in this example between 8 and 400.
#
# (In terms of Benjamini and Speed (2011), we compare W_a,l for a=margin and different l's
# We expect the best value to approximately be the median fragment length
######

# We always recommend throwing away the first few bp's, so as not
# to confound with fragmentation effects. 
margin = 5

begdata = makeGCLens(chr1$isgc,dat1$forw,sampline = sampCh1$singleLocSamp,minlen = 8,maxlens=650,margin=5,max_frag_for_loc=10)

#####
# begdata consists of:
# begdata$locs: a matrix of locations
# begdata$frag: a matrix of fragment counts
# row l in each corresponds W_{0,l} (stratification based on a GC window of size l)
#####

#####
# Computes the TV score for each row (deviation from uniformity)
####
tvs = scoreGCLens(begdata, maxlen=650, minlen =  8,scale= T)
plotGCLens(tvs,lw=2,lt=1)
best_size = which.max(tvs)-1

####
# We now use this size to generate condMean table...
# We want to make a gc line of window-2*margin
####

gcsize = best_size # I got 555 for my sample of the full chromosome
gcline = prepareGCChrom(chr1,gcsize,filename = "gcchr1")


##
# gcline[i] = isgc[i]+isgc[i+1]+....isgc[i+gcsize-1]
# We want to stratify location i based on gcline[i+margin]
# We use our sample for this, comparing the gcline[sampCh1$singleLocSamp+margin]
# to the fragment counts at these locations (already computed in sampCh1$forsamped)
##

cMeans= getCondMean(gcline[sampCh1$singleLocSamp+margin],sampCh1$forsamped,cutoff = 4,jump = 6,norm = FALSE)
# Taking jumps at jumps of 3 give smoother results.

# Similarly, we could estimate the condMean on the reverse strand
# Here, we need to do the adjustment:
# First find the 3' end (x + readlen), then go back -(margin+gcsize)
cMeans_rev= getCondMean(gcline[sampCh1$singleLocSamp+dat1$readlen-1-margin-gcsize]
  ,sampCh1$revsamped,cutoff = 4,jump = 6,norm = FALSE)
  

# And plot
par(mfcol = c(1,2))
plotCondMean(cMean = cMeans,ci = TRUE,normRange = gcsize,meanLine=TRUE,lt = 1,col=4)
plotCondMean(cMean = cMeans_rev,ci = TRUE,normRange = gcsize,meanLine=TRUE,lt = 1,col=4)

# If they seem congruent, we can add them together:
cMeans_tog = sumCondMeans(cMeans ,cMeans_rev)

# Repeat this for the other chromosomes:
# A. Create chr2, chr3, ...
# B. Create dat2, dat3, ... (data structures)
# C. Create samp2, samp3, ... (sample structures)
#    Might take sample sizes to be proportional to chromosome len.
# D. Estimate cMean2, ... based on a gcline with the
#    the same fragsize and margin. Make sure jump is the same.
# E. Add all cMean structures together using sumCondMeans


# Get predictions based on the conditional mean
# Works for each chromosome separately
forpreds = predictLine(cMeans_tog,gcline,gcsize,margin,chr1$isrep,strand = "F",readlen=dat1$readlen)
revpreds = predictLine(cMeans_tog,gcline,gcsize,margin,chr1$isrep,strand = "R",readlen=dat1$readlen)


# We can also run everything together, creating predictions based on GC window 
# (Supplying filename allows loading the gc line rather than recomputing it)
forpredsB = condMeanAndPredict(chr1,sampCh1,winstart=margin,winend=margin+gcsize,cutoff=4,gcfilename = "gcchr1",jump=6)
  

# res$pred is a single-bp rate prediction. We can pool these together to get 1K predictions
pred1K = meanLine(forpreds$preds,1000)*1000
pred1K_rev= meanLine(revpreds$preds,1000)*1000

# If we have predefined regions that are uneven, we can use
# uneven binning for the predictions
roi_begin= c(10000,300000,400000)
roi_end = c(11000,300500,400100)
roi_preds = unEvenBin(forpreds$preds,roi_begin,roi_end)

#create the 1K bincounts as well
chr1$gc1K = meanLine(chr1$isgc,1000)
chr1$map1K = meanLine(chr1$isrep,1000)
dat1$for1K = binReads(dat1$forw[,1],1000,dat1$chr_len)
dat1$rev1K = binReads(dat1$reve[,1],1000,dat1$chr_len)

# Loess is an efficient way to create predictions.
# For plotting - "uncorrLoess" smoothes the raw counts.
# For prediction - "predLoess" estimates curves on counts corrected for mappability
#                   and provides a prediction line
#
# loess1K = predLoess(dat1, chr1, binsize = "1K", workrange = 1:240000 ,span = 0.3, plot=F)

# A quick and dirty normalization based on the predicted values

norm1K = (dat1$for1K+dat1$rev1K)/(pred1K+pred1K_rev+0.01)
unnorm1K = (dat1$for1K+dat1$rev1K)/median(dat1$for1K+dat1$rev1K, na.rm=TRUE)

# Restrict our comparison to bins with high mappability (low repeats):
rang = 3000:8000
rang = rang[chr1$map1K<0.4]
par(mfcol = c(2,1))
plot(rang,norm1K[rang],cex = 0.5, pch =20)
plot(rang,unnorm1K[rang],cex = 0.5, pch =20)

# Output results
write.table(file ="chr1_for.output"
            ,data.frame(binstart=seq(1,1000000,1000),
                  binend=seq(1000,1000000,1000),
                  cntfwd = dat1$for1K[1:1000],
                  cntrev = dat1$rev1K[1:1000],
                  forpred = pred1K[1:1000],
                  revpred = pred1K_rev[1:1000],
                  repeats = chr1$map1K[1:1000],
                  gc = chr1$gc1K[1:1000]))

               
         

