
#############
#
# Maintained by Yuval Benjamini
# Adapated from Oleg Mayba, UC Berkeley
#
# This python script reads a set of reference chromosomes and generates
# all possible reads. The resulting read file can be remapped to the
# reference to identify non-unique mappings and create mappability line.
#
#############

import sys
import os
import time
import optparse

parser = optparse.OptionParser()

parser.add_option('-r','--reference_dir',
    action="store", dest="reference_dir",
    help="Name of reference dir", default="./")
parser.add_option('-o','--output_dir',
    action="store", dest="output_dir",
    help="Name of output dir", default="./")
parser.add_option('-l','--readlen',
    action="store", dest="read_length",
    help="Length of reads", default="75")
parser.add_option('-c','--chrs',
    action="store", dest="chrs",
    help="Chromosomes", default="range(1,23)+ ['X','Y','M']")

options, args = parser.parse_args()

read_length = eval(options.read_length)
chrs = eval(options.chrs)
output_dir = options.output_dir
reference_dir = options.reference_dir

def seq_from_fasta_file(fasta_file):
    fasta_seq=open(fasta_file)
    seq=''
    for line in fasta_seq:
        if line.startswith('>'):
            continue
        else:
            seq=seq+line.strip()
    fasta_seq.close()
    return seq.upper()




for chr in chrs:
    print 'chr%s' %(chr)
    chr_file=reference_dir+'chr%s.fa' %(chr)
    chr_header='>chr%s' %(chr)
    chr_seq=seq_from_fasta_file(chr_file)
    standard_bases=['A','C','G','T','N']
    nonstandard_count=0
    for i in xrange(len(chr_seq)):
        if not(chr_seq[i] in standard_bases):
           chr_seq = chr_seq[:i] + 'N' + chr_seq[i+1:] 	
           nonstandard_count+=1
    print 'Nonstandard count is %s' %(nonstandard_count)
    print 'chr length is %s' %(len(chr_seq))
    num_reads=len(chr_seq)-read_length+1
    raw_seqs_filename= output_dir + 'chr%s_seqs.raw' %(chr)
    print 'start writing seqs to file'
    raw_seqs_file=open(raw_seqs_filename,'w')
    for i in xrange(num_reads):
        raw_seq=chr_seq[i:(i+read_length)]+'\n'
        raw_seqs_file.write(raw_seq)
    raw_seqs_file.close()




