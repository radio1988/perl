# misc codes
- format conversion: fasta, fastq, gff, gtf, etc

# fasta fastq
- fastq_to_fasta.pl
- fastq_rc.pl: reverse complement


# alignment result, annotation parsing
- sam2gtf.pl: 
  use CIGAR to convert SAM to GTF, so that assembly software (such as cuffmerge) can use alignment result from long reads (such as Iso-seq, NCBI-mRNA)
  e.g.: 
    `cd sample_data` 
    `perl ../sam2gtf.pl isoseq_sam_sample.sam`
-  
