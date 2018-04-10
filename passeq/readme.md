# Code related to WTTS-seq analysis
WTTS-seq is very similar to PAS-seq 
http://www.genetics.org/content/203/2/683

# Demo Instructions:
1. `perl ./polyAsite_finder_from_sam.pl sample_data/wtts_gsnap_sample.sam`
    - find polyA site from SAM (BWA/GSNAP/GMAP produced) file
    - saved in in-house `info` format:`ref dir pos count MaxRead cigar sam-pos sam-seq`
2. `perl ./polyAsite_clustering.pl sample_data/wtts_sample.info.JGI_symbol 20`
    - find clusters of polyA sites by single linkage clustering. (if two polyA sites are close enough, they are clustered together, 20bp is the clustering distance here)
    - input: `sample_data/wtts_sample.info.JGI_symbol`, which is `info` file annotated with the help of cuffcompare and GTF/GFF genome annotation, with the help of some other scripts
    - output: 
      - `wtts_sample.info.JGI_symbol.20clust`: detailed version, format: `cluster_index scaffold direction position range count class_code longest_read_id mapping_quality cigar longest_read_length read_sequence`
      - `wtts_sample.info.JGI_symbol.20clust2`: annotated version, format: `cluster_index   scaffold        direction       position        range   count   class_code      gene_acc        gene_page       gene_symbol     gene_desc`
      - `class_code`: cuffcompare positional relationship to annotated gene, `i: intron match, =: full match, etc `,  http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/

# Full pipeline (not all scripts cleaned and uploaded yet):
* `sh pas_pipeline_from_sam.xt9.sh`
