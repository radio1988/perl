# format conversion related scripts
- format conversion: fasta, fastq, GFF, GTF, SAM, etc

## fasta fastq related
- `fastq_to_fasta.pl`: `perl fastq_to_fasta.pl sample_data/wtts_sample.fastq sample_data/wtts_sample.fasta`
- `fastq_rc.pl`: reverse complement, `perl fastq_rc.pl sample_data/wtts_sample.fastq sample_data/wtts_sample.rc.fastq`


## alignment result, annotation parsing related
- `sam2gtf.pl`
	- `perl sam2gtf.pl sample_data/isoseq_sam_sample.sam` 
	- use CIGAR to convert SAM to GTF, so that assembly software (such as cuffmerge) can use alignment result from long reads (such as Iso-seq, NCBI-mRNA)

- `gff3_to_gtf.pl`
	- `perl gff3_to_gtf.pl sample_data/Xentr7_2_Stable.mRNA.sorted.scaffold_5389.gff3` 
