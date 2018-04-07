#!/usr/bin/perl
#Rui Li, Animal Sciences Department, Washington State University;
#PhD student of Dr. Zhihua Jiang
#Email: liruiradiant+perl@gmail.com

use strict; use warnings;
use READ_GFF;
print"
input: Xenbase7.2.gff3
output: gtf format for gtf_splicesites in gsnap

usage:<pl> <*.gff3>
2014-11-21
";

die "cmd err\n" unless @ARGV == 1;
open(IN,$ARGV[0]) or die "err reading\n";
open(OUT,">$ARGV[0].gtf") or die "err outputing\n";

my@in = (<IN>);

my%pacid_2_geneid;
my%pacid_2_transcriptid;

foreach(@in){
	if(/\tmRNA\t/){
		my%gff = READ_GFF::XENBASE($_);
		$pacid_2_geneid{ $gff{pacid} } = $gff{parent} ;
		$pacid_2_transcriptid{ $gff{pacid}} = $gff{name};
	}
}

#foreach(keys %pacid_2_geneid){	print"gene:$_ - $pacid_2_geneid{$_}\n";next}
#foreach(keys %pacid_2_transcriptid){	print"mRNA:$_ - $pacid_2_transcriptid{$_}\n";next}

foreach(@in){
	chomp;
	if(/^#/){
		print OUT "$_\n";
		next;
	}#for header of gff

	my%gff = READ_GFF::XENBASE($_);
	
	my@ele = split /\t/,$_;
	my$unchanged = join "\t", @ele[0..7];

	if(/\tgene\t/){
		print OUT "$unchanged\tgene_id \"$gff{name}\";\n";		
		next;
	}
	
	my$geneid = $pacid_2_geneid{ $gff{pacid} } ;
	my$transcriptid = $pacid_2_transcriptid{ $gff{pacid} } ;
	print OUT "$unchanged\tgene_id \"$geneid\";transcript_id \"$transcriptid\";\n";
}
