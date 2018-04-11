#!/usr/bin/perl
# Rui Li, Animal Sciences Department, Washington State University;
# PhD student of Dr. Zhihua Jiang
# Email: liruiradiant+perl@gmail.com

use strict; 
use READ_SAM;

print"
2015-08-06
sam-position + cigar => gff
usage:<pl> <*.sam>
output: *.sam.gtf

Note: 
you can ignore warning of `NM missing` for CLC-SAM file, which is no-standard SAM
The output for CLC-SAM will still be correct
";

die "usage: <sam2gtf.pl> <sam>\n" unless @ARGV == 1;

open(IN,$ARGV[0]) or die;
open(OUT,">$ARGV[0].gtf") or die;


while(<IN>){
	my%sam = READ_SAM::General($_);
	if($sam{flag} == 4){next}
	my$direction = "NA";
	if ($sam{flag} & 16){$direction = "-"}else{$direction = "+"}
	my$gff = Cigar_2_Gff($sam{rname},$sam{pos},$direction,$sam{cigar},$sam{qname}); 
}

sub Cigar_2_Gff{
	my($scaffold,$position,$direction,$cigar,$id)= @_;
#       print"$position\t$cigar\t$direction\t$id\t$scaffold\n";
#       print"$cigar\n";
	my@segments = split /N/,$cigar;
	my$exon_number;
	my$exon_start=$position;
	foreach(@segments){ #part of the cigar
		$exon_number++;
		my($length) = Cigar_2_Segment($_);
		my$intron;
		if($_ =~ /\d+$/){($intron) =  $_ =~ /(\d+$)/}else{$intron=0};
		print OUT "$scaffold\tsam\texon\t$exon_start\t",$exon_start+$length-1,"\t.\t$direction\t.\tgene_id \"$id\"; transcript_id \"$id\"; exon_number \"$exon_number\";\n";
		$exon_start += ($length-1+$intron+1);
	}
}


sub Cigar_2_Segment{
        my%sam;
        $sam{cigar} = shift;
        my@cigar_S = $sam{cigar} =~ /(\d+)S/g;
        my@cigar_H = $sam{cigar} =~ /(\d+)H/g;
        my@cigar_M = $sam{cigar} =~ /(\d+)M/g;
        my@cigar_I = $sam{cigar} =~ /(\d+)I/g;
        my@cigar_D = $sam{cigar} =~ /(\d+)D/g;
        my@cigar_N = $sam{cigar} =~ /(\d+)N/g;
        my@cigar_X = $sam{cigar} =~ /(\d+)X/g;
        my@cigar_euqal = $sam{cigar} =~ /(\d+)=/g;
        my($mcount,$icount,$dcount,$scount,$hcount,$ncount,$Xcount,$Equal_count) = (0,0,0,0,0,0,0,0);
        foreach(@cigar_S){$scount += $_}
        foreach(@cigar_H){$hcount += $_}
        foreach(@cigar_M){$mcount += $_}
        foreach(@cigar_I){$icount += $_}
        foreach(@cigar_D){$dcount += $_}
        foreach(@cigar_N){$ncount += $_}
        foreach(@cigar_X){$Xcount += $_}
        foreach(@cigar_euqal){$Equal_count += $_}
        $sam{cigar_length} = $mcount + $dcount + $Xcount + $Equal_count;
        return($sam{cigar_length});


}