# format conversion related scripts
- format conversion: fasta, fastq, GFF, GTF, SAM, etc

## fasta fastq related
- `fastq_to_fasta.pl`
- `fastq_rc.pl`: reverse complement


## alignment result, annotation parsing related
- `sam2gtf.pl`: 
  - use CIGAR to convert SAM to GTF, 
  - so that assembly software (such as cuffmerge) can use alignment result from long reads (such as Iso-seq, NCBI-mRNA)
  - e.g.: 
    1. `cd sample_data` 
    2. `perl ../sam2gtf.pl isoseq_sam_sample.sam`

- `` 
