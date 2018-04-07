#!/usr/bin/perl
#Rui Li, Animal Sciences Department, Washington State University;
#PhD student of Dr. Zhihua Jiang
#Email: liruiradiant+perl@gmail.com

use strict; use warnings;
use READ_GFF;
print"
input: #!gff-spec-version 1.21
output: ensembl gtf format for snpEFF

usage:<gff1.21_to_gtf.pl> <xxx.gff>
2016-01-22
";

die "cmd err\n" unless @ARGV == 1;
open(IN,$ARGV[0]) or die "err reading\n";
open(OUT,">$ARGV[0].gtf") or die "err outputing\n";

my@in = (<IN>);
close IN;

my%transid2geneid;
my%pacid_2_transcript_id;

foreach(@in){
	if(/Parent=/){
		my%gff = READ_GFF::XENBASE($_);
		$transid2geneid{ $gff{ID} } = $gff{parent} ;
		$pacid_2_transcript_id{ $gff{ID}} = $gff{ID};
	}
}

#foreach(keys %transid2geneid){	print"gene:$_ - $transid2geneid{$_}\n";next}
#foreach(keys %pacid_2_transcript_id){	print"mRNA:$_ - $pacid_2_transcript_id{$_}\n";next}

foreach(@in){
	chomp;

	if(/^#/){
		#print OUT "$_\n";
		next;
	}#for header of gff

	my%gff = READ_GFF::GFF121($_);
	my@ele = split /\t/,$_;
	my$unchanged = join "\t", @ele[0..7];

	#for tRNA
	if(/\ttRNAscan-SE\t/) {next; }




#for mRNA
#	if(/\t(mRNA|transcript)\t/) {
#		my$geneid = $gff{parent};
#		my$transcript_id = $gff{ID};
#		my$transcript_name = $gff{name};
#		my$gene_name = $gff{gene};
#		print OUT "$unchanged\tgene_id \"$geneid\";transcript_id \"$transcript_id\";gene_name \"$gene_name\";transcript_name \"$transcript_name\";\n";
#		next;
#	}

	if(/\texon\t/){
		my$transcript_id = $gff{parent};
		my$geneid = $transid2geneid{$transcript_id};
		my$transcript_name = $gff{transcript_id};
		my$gene_name = $gff{gene};
		print OUT "$unchanged\tgene_id \"$geneid\"; transcript_id \"$transcript_id\"; gene_name \"$gene_name\"; transcript_name \"$transcript_name\";\n";
		next;
	}
	
	if(/\tCDS\t/){
		my$transcript_id = $gff{parent};
		my$geneid = $transid2geneid{$transcript_id};
		my$protein_name = $gff{name};
		my$gene_name = $gff{gene};
		print OUT "$unchanged\tgene_id \"$geneid\"; transcript_id \"$transcript_id\"; gene_name \"$gene_name\"; protein_name \"$protein_name\";\n";
		next;
	}
	
	
	
}
