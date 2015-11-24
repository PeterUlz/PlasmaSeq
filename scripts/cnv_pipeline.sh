#! /bin/bash
#
# CNV calling Pipeline (see Baslan et al. (Nat Prot., 2012)
#

project=$1
fastq=$2
proj_dir=$PWD/$project
CGHweb_dir=$proj_dir/CGHResults
mkdir $proj_dir
echo $fastq `date` `md5sum $fastq`> $proj_dir/${project}.data
#align and output sam file
/home/peter/Software/Alignment/bwa-0.6.2/bwa aln -f ${1}.aln -t 20 ~/RefSeq/Navin_CNV/hg19/hg19.navin $2
/home/peter/Software/Alignment/bwa-0.6.2/bwa samse -f ${1}.sam ~/RefSeq/Navin_CNV/hg19/hg19.navin ${1}.aln $2 

#create bam file, sort and remove duplicates
samtools view -S -b -o ${1}.bam ${1}.sam
rm ${1}.sam
rm ${1}.aln
samtools rmdup -s ${1}.bam ${1}.rmdup.bam
samtools sort ${1}.rmdup.bam ${1}.rmdup.sorted
samtools view ${1}.rmdup.sorted.bam > $proj_dir/${1}.sam
rm ${1}.bam
rm ${1}.rmdup.bam
rm ${1}.rmdup.sorted.bam

#count reads in bins
/home/peter/Pipeline/CNV_calling/Breast/count_reads_in_bins_cnv.py $proj_dir/${1}.sam $proj_dir/${1}.bincounts $proj_dir/${1}.stats
cd $proj_dir

rm ${1}.sam

#normalize GC create plots and segment
echo "cbs.segment01(indir=\".\", outdir=\".\", bad.bins=\"/home/peter/RefSeq/Navin_CNV/hg19.50k.k50.bad.bins.txt\", varbin.gc=\"/home/peter/RefSeq/Navin_CNV/hg19.new_sorted.gc_count.txt\", varbin.data=\"${1}.bincounts\", sample.name=\"$1\", alt.sample.name=\"\", alpha=0.05, nperm=1000, undo.SD=1.0, min.width=5)" > /home/peter/Pipeline/CNV_calling/Breast/tmp.R

cat /home/peter/Pipeline/CNV_calling/Breast/CNV_postprocessing.cghweb.r /home/peter/Pipeline/CNV_calling/Breast/tmp.R > /home/peter/Pipeline/CNV_calling/Breast/execute.R

R CMD BATCH /home/peter/Pipeline/CNV_calling/Breast/execute.R ${1}.out

rm /home/peter/Pipeline/CNV_calling/Breast/execute.R
rm /home/peter/Pipeline/CNV_calling/Breast/tmp.R

gunzip $CGHweb_dir/Table_of_aCGH_smoothed_profiles.txt.gz
cat $CGHweb_dir/Table_of_aCGH_smoothed_profiles.txt > $1.segmented

/home/peter/Pipeline/CNV_calling/Breast/get_segments.pl $1.segmented $1.segments

cat /home/peter/Pipeline/CNV_calling/Breast/makeGraph.R | R --slave --args 2 $1.segmented
cat /home/peter/Pipeline/CNV_calling/Breast/makeGraph_each_chr.R | R --slave --args 2 $1.segmented
cat /home/peter/Pipeline/CNV_calling/Breast/makeGraph_segments.R | R --slave --args 2 $1.segmented
cat /home/peter/Pipeline/CNV_calling/Breast/makeGraph_segments_each_chr.R | R --slave --args 2 $1.segmented
cat /home/peter/Pipeline/CNV_calling/Breast/makeGraph_linear_chrlen_log2.R | R --slave --args 2 $1.segmented
cat /home/peter/Pipeline/CNV_calling/Breast/makeGraph_linear_chrlen_segments.R | R --slave --args 2 $1.segmented
cat /home/peter/Pipeline/CNV_calling/Breast/makeGraph_linear_chrlen_blue.R | R --slave --args 2 $1.segments

/home/peter/Pipeline/CNV_calling/Breast/segmental_GC_z-scores.pl $1.segments ${1}.corrected.bincounts ${1}.segmented.zscores.txt

#Focal Amplification mapping
mkdir FocalAmps
cat /home/peter/Pipeline/CNV_calling/Breast/FocalAmplifications_fromPlasmaSeq.R  | R --slave --args $proj_dir/$1.segmented $proj_dir/$1.segments $proj_dir/FocalAmps/${1} Breast 100
