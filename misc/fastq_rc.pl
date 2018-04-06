#!/usr/bin/perl
#Rui Li, Animal Sciences Department, Washington State University;
#PhD student of Dr. Zhihua Jiang
#Email: liruiradiant+perl@gmail.com

use strict; use warnings;
print'
usage: <fastq_to_fasta.pl> <inname.fastq> <outname.fastq>
example: perl fastq_to_fasta.pl test.fastq test.rc.fastq
reverse compliment fastq file
';
die "err in cmd\n" unless @ARGV == 2;
open (IN, "<$ARGV[0]") or die "err reading input\n";
open (OUT, ">$ARGV[1]") or die "err writing output\n";

my%fq;my$i = 0;my $fq_count;
while(<IN>){
	$_ =~ s/\r//;
	$_ =~ s/\n//;
	if(/^@/){$i = 0}
	$i ++;

	if($i == 1) {$fq{name} = $_}
	if($i == 2) {$fq{seq} =  $_}
	if($i == 4) {$fq{quali} = $_}

	if($i == 4){
		$fq{seq} =~ tr/ATGCatgc/TACGtacg/;
		$fq{seq} = reverse $fq{seq};
		$fq{quali} = reverse $fq{quali};
		print OUT "$fq{name}\n$fq{seq}\n+\n$fq{quali}\n";
		%fq = ();
		$fq_count++;
	}
}

print "there are $fq_count seqs in total\n\n"
