# PlasmaSeq
Analyze CNVs of ctDNA

*Attention: Large hg19 reference index files (with pseudo-autosomal region masked) for BWA are not in this commit*

CNVs of ctDNA low-coverage WGS sequencing are analyzed in several steps

 1) Merge FastQ Files for different lanes (NextSeq only)  
 2) Align FastQ to modified hg19 (pseudo-autosomal region of chrY masked)  
 3) Remove PCR duplicates  
 4) Count reads in (50,000) predefined genomic bins (average length ~56kbp, each bin containing equal amount of mappable positions)  
 5) Normalize by mean read count  
 6) Normalize by GC-content (Lowess-smoothing)  
 7) Normalize by mean (normalized) read count per bin of healthy controls  
 8) Segment using combination of GLAD and CBS (provided by CGHWeb package in R)  
 9) Create plots  
10) Identify Focal amplifications and deletions
