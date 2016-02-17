#! /usr/bin/python

# Pipeline for CNV analysis of Plasma-Seq data

#version 0.3: run R in isolation mode
#version 0.2: incorporate nextseq support
from subprocess import call
import sys
import argparse
import os
import time
import shutil


script_dir = os.path.dirname(os.path.realpath(__file__))
# Parse command line arguments ###################################################################################
parser = argparse.ArgumentParser(description='Analyze read depth in comparison to transcription end')
parser.add_argument('-fq','--fastq', dest='fastq_file', 
                   help='Fastq file (for NextSeq: Specify Lane 1 Fastq File)',required=True)
parser.add_argument('-s','--sample-name', dest='name',
                   help='Sample name to be used subsequently',required=True)
parser.add_argument('-g','--gender', dest='gender',
                   help='Gender of the sample [either m or f]',required=True, choices=["m","f"])
parser.add_argument('-o','--out-dir', dest='outdir',
                   help='Output Directory [default .]',default=".")
parser.add_argument('-k','--keep-temp', dest='keep',
                   help='Keep temporary files',action="store_true")
parser.add_argument('-m','--machine', dest='machine',
                   help='Sequencing machine [miseq|nextseq]',required=True, choices=["miseq","nextseq"])
parser.add_argument('-skipmerge','--skip-merge', dest='skip_merge',
                   help='Skip merging of FastQ files for NextSeq data',action="store_true")
parser.add_argument('-t','--threads', dest='threads',
                   help='No. threads for alignment [default: 1]',type=int,default=1)

args = parser.parse_args()
print time.strftime("%d/%m/%Y:  %H:%M:%S  : Starting Analysis")
print "  Fastq file: ",args.fastq_file
print "  Sample name file: ",args.name
print "  Output Directory: ",args.outdir
print "  Machine: ",args.machine
print "  Gender: ",args.gender
if args.keep:
   print "  Keeping temporary files"

proj_dir = args.outdir+"/"+args.name

########################################################################################################
# Step 0.5 merge FastQ Files if nextseq is used
if args.machine == "nextseq" and not args.skip_merge:
    if "L001_R1_001" not in args.fastq_file:
        print "Please specify FastQ File of Lane1 for NextSeq data"
        sys.exit(1)
    lane1_file = args.fastq_file
    lane2_file = args.fastq_file[:-17]+"2_R1_001.fastq.gz"
    lane3_file = args.fastq_file[:-17]+"3_R1_001.fastq.gz"
    lane4_file = args.fastq_file[:-17]+"4_R1_001.fastq.gz"
    outfilename = args.fastq_file[:-21]+".fastq.gz"
    with open(outfilename, 'wb') as outfile:
        for filename in [lane1_file,lane2_file,lane3_file,lane4_file]:
            with open(filename, 'rb') as readfile:
                shutil.copyfileobj(readfile, outfile)
    args.fastq_file = outfilename
########################################################################################################
# Step1 create directory and create MD5 file of input
return_val=call(["mkdir",proj_dir])
if return_val != 0:
    print "Cannot create directory "+proj_dir 
    sys.exit(1)

ERR_LOG = open(proj_dir+"/"+args.name+".err.log","w")
OUT_LOG = open(proj_dir+"/"+args.name+".out.log","w")

OUT_LOG.write(time.strftime("%d/%m/%Y:  %H:%M:%S  : Starting Analysis"))
OUT_LOG.write("Fastq file: "+args.fastq_file+"\n")
OUT_LOG.write("Sample name file: "+args.name+"\n")
OUT_LOG.write("Output Directory: "+args.outdir+"\n")
OUT_LOG.write("Machine: "+args.machine+"\n")
OUT_LOG.write("Gender: "+args.gender+"\n")
if args.keep:
   OUT_LOG.write("Keeping temporary files"+"\n")
OUT_LOG.write(time.strftime("\n%d/%m/%Y:  %H:%M:%S  : Step1 create directory and create MD5 file of input\n"))
print time.strftime("%d/%m/%Y:  %H:%M:%S  : Step1 create directory and create MD5 file of input")

MD5_OUTPUT = open(proj_dir +"/"+args.name+".md5","w")
call(["md5sum",args.fastq_file],stdout=MD5_OUTPUT,stderr=ERR_LOG)
MD5_OUTPUT.close()

########################################################################################################
# Step2 alignment
OUT_LOG.write(time.strftime("\n%d/%m/%Y:  %H:%M:%S  : Step2 alignment\n"))
print time.strftime("%d/%m/%Y:  %H:%M:%S  : Step2 alignment")
aln_value=call([script_dir+"/software/bwa","aln","-f",proj_dir+"/"+args.name+".aln","-t",str(args.threads),script_dir+"/ref/hg19.navin",args.fastq_file],stderr=ERR_LOG,stdout=OUT_LOG)
if aln_value != 0:
    print "Error in alignment" 
    sys.exit()
aln_value=call([script_dir+"/software/bwa","samse","-f",proj_dir+"/"+args.name+".sam",script_dir+"/ref/hg19.navin",proj_dir+"/"+args.name+".aln",args.fastq_file],stderr=ERR_LOG,stdout=OUT_LOG)
if aln_value != 0:
    print "Error in alignment" 
    sys.exit()

######################################################################################################## 
# Step3 create bam file, sort and remove duplicates
OUT_LOG.write(time.strftime("\n%d/%m/%Y:  %H:%M:%S  : Step3 create bam file, sort and remove duplicates\n"))
print time.strftime("%d/%m/%Y:  %H:%M:%S  : Step3 create bam file, sort and remove duplicates")
call([script_dir+"/software/samtools","view","-S","-b","-o",proj_dir+"/"+args.name+".bam",proj_dir+"/"+args.name+".sam"],stderr=ERR_LOG,stdout=OUT_LOG)
if not args.keep:
    call(["rm",proj_dir+"/"+args.name+".sam"])
    call(["rm",proj_dir+"/"+args.name+".aln"])
call([script_dir+"/software/samtools","rmdup","-s",proj_dir+"/"+args.name+".bam",proj_dir+"/"+args.name+".rmdup.bam"],stderr=ERR_LOG,stdout=OUT_LOG)
call([script_dir+"/software/samtools","sort",proj_dir+"/"+args.name+".rmdup.bam",proj_dir+"/"+args.name+".rmdup.sorted"],stderr=ERR_LOG,stdout=OUT_LOG)
call([script_dir+"/software/samtools","view","-o",proj_dir+"/"+args.name+".rmdup.sorted.sam",proj_dir+"/"+args.name+".rmdup.sorted.bam"],stderr=ERR_LOG,stdout=OUT_LOG)

if not args.keep:
    call(["rm",proj_dir+"/"+args.name+".bam"],stderr=ERR_LOG,stdout=OUT_LOG)
    call(["rm",proj_dir+"/"+args.name+".rmdup.bam"],stderr=ERR_LOG,stdout=OUT_LOG)
    call(["rm",proj_dir+"/"+args.name+".rmdup.sorted.bam"],stderr=ERR_LOG,stdout=OUT_LOG)


######################################################################################################## 
# Step4 count reads in predefined genomic bins
OUT_LOG.write(time.strftime("\n%d/%m/%Y:  %H:%M:%S  : Step4 count reads in predefined genomic bins\n"))
print time.strftime("%d/%m/%Y:  %H:%M:%S  : Step4 count reads in predefined genomic bins")
call([script_dir+"/scripts/count_reads_in_bins_cnv.py",proj_dir+"/"+args.name+".rmdup.sorted.sam",proj_dir+"/"+args.name+".bincounts",proj_dir+"/"+args.name+".stats",script_dir+"/ref"],stderr=ERR_LOG,stdout=OUT_LOG)

if not args.keep:
    call(["rm",proj_dir+"/"+args.name+".rmdup.sorted.sam"],stderr=ERR_LOG,stdout=OUT_LOG)

######################################################################################################## 
# Step5 Normalization in R (GC-correction, normalization with means of controls, etc...) and segmentation
OUT_LOG.write(time.strftime("\n%d/%m/%Y:  %H:%M:%S  : Step5 Normalization in R (GC-correction, normalization with means of controls, etc...) and segmentation\n"))
print time.strftime("%d/%m/%Y:  %H:%M:%S  : Step5 Normalization in R (GC-correction, normalization with means of controls, etc...) and segmentation")
TMP_R = open(proj_dir+"/"+args.name+".tmp.R","w")
if args.machine == "miseq":
    if args.gender == "m":
        TMP_R.write("cbs.segment01(indir=\".\", outdir=\""+proj_dir+"/CGHResults\", bad.bins=\""+script_dir+"/ref/hg19.50k.k50.bad.bins.txt\","+
        "varbin.gc=\""+script_dir+"/ref/hg19.new_sorted.gc_count.txt\", varbin.data=\""+proj_dir+"/"+args.name+".bincounts"+"\", sample.name=\""+args.name+"\", "+
        "alt.sample.name=\"\", alpha=0.05, nperm=1000, undo.SD=1.0, min.width=5,controls_file=\""+script_dir+"/ref/Kontrollen_male.bincount.txt\",sample.dir=\""+proj_dir+"\")\n")
    else:
        TMP_R.write("cbs.segment01(indir=\".\", outdir=\""+proj_dir+"/CGHResults\", bad.bins=\""+script_dir+"/ref/hg19.50k.k50.bad.bins.txt\","+
        "varbin.gc=\""+script_dir+"/ref/hg19.new_sorted.gc_count.txt\", varbin.data=\""+proj_dir+"/"+args.name+".bincounts"+"\", sample.name=\""+args.name+"\", "+
        "alt.sample.name=\"\", alpha=0.05, nperm=1000, undo.SD=1.0, min.width=5,controls_file=\""+script_dir+"/ref/Kontrollen_female.bincount.txt\",sample.dir=\""+proj_dir+"\")\n")
if args.machine == "nextseq":
    if args.gender == "m":
        TMP_R.write("cbs.segment01(indir=\".\", outdir=\""+proj_dir+"/CGHResults\", bad.bins=\""+script_dir+"/ref/hg19.50k.k50.bad.bins.txt\","+
        "varbin.gc=\""+script_dir+"/ref/hg19.new_sorted.gc_count.txt\", varbin.data=\""+proj_dir+"/"+args.name+".bincounts"+"\", sample.name=\""+args.name+"\", "+
        "alt.sample.name=\"\", alpha=0.05, nperm=1000, undo.SD=1.0, min.width=5,controls_file=\""+script_dir+"/ref/Kontrollen_male.bincount_nextseq.txt\",sample.dir=\""+proj_dir+"\")\n")
    else:
    # TODO: edit to include nextseq female controls
        TMP_R.write("cbs.segment01(indir=\".\", outdir=\""+proj_dir+"/CGHResults\", bad.bins=\""+script_dir+"/ref/hg19.50k.k50.bad.bins.txt\","+
        "varbin.gc=\""+script_dir+"/ref/hg19.new_sorted.gc_count.txt\", varbin.data=\""+proj_dir+"/"+args.name+".bincounts"+"\", sample.name=\""+args.name+"\", "+
        "alt.sample.name=\"\", alpha=0.05, nperm=1000, undo.SD=1.0, min.width=5,controls_file=\""+script_dir+"/ref/Kontrollen_female.bincount.txt\",sample.dir=\""+proj_dir+"\")\n")
TMP_R.close()

R_SCRIPT = open(proj_dir+"/"+args.name+".script.R","w")
COMMON_R = open(script_dir+"/scripts/CNV_postprocessing.cghweb.r","r")
for line in COMMON_R.readlines():
    R_SCRIPT.write(line)
COMMON_R.close()
TMP_R = open(proj_dir+"/"+args.name+".tmp.R","r")
for line in TMP_R.readlines():
    R_SCRIPT.write(line)
TMP_R.close()
R_SCRIPT.close()
os.environ["R_LIBS_USER"]=script_dir+"/R/x86_64-pc-linux-gnu-library/2.14"
r_loc=script_dir+"/R/R"
call([r_loc,"CMD","BATCH",proj_dir+"/"+args.name+".script.R",proj_dir+"/"+args.name+".script.R.out"],stderr=ERR_LOG,stdout=OUT_LOG)
call(["gunzip",proj_dir+"/CGHResults/Table_of_aCGH_smoothed_profiles.txt.gz"],stderr=ERR_LOG,stdout=OUT_LOG)
call(["mv",proj_dir+"/CGHResults/Table_of_aCGH_smoothed_profiles.txt",proj_dir+"/"+args.name+".segmented"],stderr=ERR_LOG,stdout=OUT_LOG)
call([script_dir+"/scripts/get_segments.pl",proj_dir+"/"+args.name+".segmented",proj_dir+"/"+args.name+".segments"],stderr=ERR_LOG,stdout=OUT_LOG)

######################################################################################################## 
# Step6 Create Plots in R
OUT_LOG.write(time.strftime("\n%d/%m/%Y:  %H:%M:%S  : Step6 Create Plots in R\n"))
print time.strftime("%d/%m/%Y:  %H:%M:%S  : Step6 Create Plots in R")
INFILE1 = open(script_dir+"/scripts/makeGraph.R","r")
call([r_loc,"--slave","--args",proj_dir+"/"+args.name+".segmented"],stdin=INFILE1,stderr=ERR_LOG,stdout=OUT_LOG)
INFILE1.close()

INFILE2 = open(script_dir+"/scripts/makeGraph_each_chr.R","r")
call([r_loc,"--slave","--args",proj_dir+"/"+args.name+".segmented"],stdin=INFILE2,stderr=ERR_LOG,stdout=OUT_LOG)
INFILE2.close()

INFILE3 = open(script_dir+"/scripts/makeGraph_segments.R","r")
call([r_loc,"--slave","--args",proj_dir+"/"+args.name+".segmented"],stdin=INFILE3,stderr=ERR_LOG,stdout=OUT_LOG)
INFILE3.close()

INFILE4 = open(script_dir+"/scripts/makeGraph_segments_each_chr.R","r")
call([r_loc,"--slave","--args",proj_dir+"/"+args.name+".segmented"],stdin=INFILE4,stderr=ERR_LOG,stdout=OUT_LOG)
INFILE4.close()

INFILE5 = open(script_dir+"/scripts/makeGraph_linear_chrlen_log2.R","r")
call([r_loc,"--slave","--args",proj_dir+"/"+args.name+".segmented"],stdin=INFILE5,stderr=ERR_LOG,stdout=OUT_LOG)
INFILE5.close()

INFILE6 = open(script_dir+"/scripts/makeGraph_linear_chrlen_segments.R","r")
call([r_loc,"--slave","--args",proj_dir+"/"+args.name+".segmented"],stdin=INFILE6,stderr=ERR_LOG,stdout=OUT_LOG)
INFILE6.close()

INFILE7 = open(script_dir+"/scripts/makeGraph_linear_chrlen_blue.R","r")
call([r_loc,"--slave","--args",proj_dir+"/"+args.name+".segments"],stdin=INFILE7,stderr=ERR_LOG,stdout=OUT_LOG)
INFILE7.close()


######################################################################################################## 
# Step7 Calculate Z-scores
OUT_LOG.write(time.strftime("\n%d/%m/%Y:  %H:%M:%S  : Step7 Calculate Z-scores\n"))
print time.strftime("%d/%m/%Y:  %H:%M:%S  : Step7 Calculate Z-scores")
if args.machine == "miseq":
    if args.gender=="m":
        call([script_dir+"/scripts/segmental_GC_z-scores.pl",proj_dir+"/"+args.name+".segments",proj_dir+"/"+args.name+".corrected.bincounts",proj_dir+"/"+args.name+".segmented.zscores.txt",script_dir+"/ref/GC_corrected_bincounts_male"],stderr=ERR_LOG,stdout=OUT_LOG)
    else:
        call([script_dir+"/scripts/segmental_GC_z-scores.pl",proj_dir+"/"+args.name+".segments",proj_dir+"/"+args.name+".corrected.bincounts",proj_dir+"/"+args.name+".segmented.zscores.txt",script_dir+"/ref/GC_corrected_bincounts_female"],stderr=ERR_LOG,stdout=OUT_LOG)
if args.machine == "nextseq":
    if args.gender=="m":
        call([script_dir+"/scripts/segmental_GC_z-scores.pl",proj_dir+"/"+args.name+".segments",proj_dir+"/"+args.name+".corrected.bincounts",proj_dir+"/"+args.name+".segmented.zscores.txt",script_dir+"/ref/GC_corrected_bincounts_male_nextseq"],stderr=ERR_LOG,stdout=OUT_LOG)
    else:
    #TODO: update command when female controls are finished analyzing
        call([script_dir+"/scripts/segmental_GC_z-scores.pl",proj_dir+"/"+args.name+".segments",proj_dir+"/"+args.name+".corrected.bincounts",proj_dir+"/"+args.name+".segmented.zscores.txt",script_dir+"/ref/GC_corrected_bincounts_female"],stderr=ERR_LOG,stdout=OUT_LOG)


######################################################################################################## 
# Step8 Call Focal Amplifications
OUT_LOG.write(time.strftime("\n%d/%m/%Y:  %H:%M:%S  : Step8 Call Focal Amplifications\n"))
print time.strftime("%d/%m/%Y:  %H:%M:%S  : Step8 Call Focal Amplifications")
call(["mkdir",proj_dir+"/FocalAmps"],stderr=ERR_LOG,stdout=OUT_LOG)
FOCAL_SCRIPT=open(script_dir+"/scripts/FocalAmplifications_fromPlasmaSeq.R","r")
if args.gender == "m":
    call(["R","--slave","--args",proj_dir+"/"+args.name+".segmented",proj_dir+"/"+args.name+".segments",proj_dir+"/FocalAmps/"+args.name,"Prostate","100",script_dir+"/scripts/Ref"],stdin=FOCAL_SCRIPT,stderr=ERR_LOG,stdout=OUT_LOG)
else:
    call(["R","--slave","--args",proj_dir+"/"+args.name+".segmented",proj_dir+"/"+args.name+".segments",proj_dir+"/FocalAmps/"+args.name,"Breast","100",script_dir+"/scripts/Ref"],stdin=FOCAL_SCRIPT,stderr=ERR_LOG,stdout=OUT_LOG)
FOCAL_SCRIPT.close()

######################################################################################################## 
OUT_LOG.write(time.strftime("\n%d/%m/%Y:  %H:%M:%S  : Finished Analysis\n"))
print time.strftime("%d/%m/%Y:  %H:%M:%S  : Finished Analysis")
ERR_LOG.close()
OUT_LOG.close()

