#!/usr/bin/perl
#Rui Li, Animal Sciences Department, Washington State University;
#PhD student of Dr. Zhihua Jiang
#Email: liruiradiant+perl@gmail.com

use strict; use warnings;
print'
usage: <pl> <fq>
reverse compliment fastq file
';
die "err in cmd\n" unless @ARGV == 1;
open (OUT, ">$ARGV[0].rc") or die "err output\n";

my%fq;my$i = 0;my $fq_count;
while(<>){
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
