# code related to WTTS-seq analysis
WTTS-seq is very similar to PAS-seq 
http://www.genetics.org/content/203/2/683

# instructions
1. `polyAsite_finder_from_sam.pl`: find polyA site from SAM (BWA/GSNAP/GMAP produced) file and save it in in-house `info` format
   - usage: `perl polyAsite_finder_from_sam.pl sample_data/wtts_gsnap_sample.sam`
   - if the alignment don't match criteria, polyA site will not be called
   - `xxx.PAS.sam`: alignement contributed to polyA sites
   - `xxx.NoPAS.sam`: alignment not matching calling polyA site criteria 
   - `xxx.sam.PAS.info`: the output file, containing information about polyA sites
     - (ref dir pos count MaxRead cigar sam-pos sam-seq)

2. `polyAsite_clustering.pl`: find clusters of polyA sites by single linkage clustering (if two polyA sites are close enough, they are clustered together)
   - usage: `perl polyAsite_clustering.pl sample_data/wtts_sample.info.JGI_symbol 20`  
     - 20bp is the example clustering distance 
   - input: `sample_data/wtts_sample.info.JGI_symbol`, `info` file annotated with the help of cuffcompare and GTF/GFF genome annotation with the help of some other scripts
   - output: 
     - `wtts_sample.info.JGI_symbol.20clust`: detailed version 
       - (cluster_index scaffold direction position range count class_code longest_read_id mapping_quality cigar longest_read_length read_sequence)
     - `wtts_sample.info.JGI_symbol.20clust2`: annotated version 
       - (cluster_index   scaffold        direction       position        range   count   class_code      gene_acc        gene_page       gene_symbol     gene_desc)
       - `class_code`: positional relationship to annotated gene

