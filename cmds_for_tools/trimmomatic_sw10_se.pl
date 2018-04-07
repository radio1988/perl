#!/usr/bin/perl
use strict; use warnings;
print "usage: <trimmomatic_sw10.pl> <suffix>\n";
die"usage: <trimmomatic_sw10.pl> <suffix>" unless @ARGV == 1;
`for file in *$ARGV[0];do echo \$file;java -jar ~/bin-from032014/Trimmomatic-0.33/trimmomatic-0.33.jar SE -threads 2 \$file \$file.trimmo10 SLIDINGWINDOW:4:10 MINLEN:16 TOPHRED33;done`;
`rename 's/fq.trimmo10/trimmo10.fq/' *.trimmo10`;
