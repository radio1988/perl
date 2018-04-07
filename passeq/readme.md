# code related to WTTS-seq analysis
WTTS-seq is very similar to PAS-seq 
http://www.genetics.org/content/203/2/683

# instructions
1. `polyAsite_finder_from_sam.pl`: find polyA site from SAM (GSNAP produced) file and save it in `info` format
  - if the alignment don't match criteria, polyA site will not be called
  - `xxx.PAS.sam`: alignement contributed to polyA sites
  - `xxx.NoPAS.sam`: alignment not matching calling polyA site criteria 
2. ``

