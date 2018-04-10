#!/usr/bin/perl
# Rui Li, Animal Sciences Department, Washington State University;
# PhD student of Dr. Zhihua Jiang
# Email: liruiradiant+perl@gmail.com

# Input: 
#	sam file of aligned WTTS-seq reads 
# 	designed for bwa-sam, gmap-sam, gsnap-sam, other SAM files may also work

# Output: 
# 	xx.info:	file containing polyA_site infomation, which is used by downstream analysis
# 	xx.PAS.sam:   sam entries with PAS called
# 	xx.NoPAS.sam: sam entries failed to contribute to PAS calling

# PAS calling:
#   The first aligned base in the input is considered the polyA-site: 
#   Read format: TTTTTXXXXXXXXXXXXXXX assumed
#   keep the first alignment if multiple alignment of the same read exists
#   multiple PAS for one read will be recorded in *.2+ file
#   -/+ refer to the direction of gsnap mapping output
#   10bp is the max overhang of mapping

# Log:
#   2014/12/17 set mini-length for PAS-Exon
#   2014/12/18 filter simple sequence for PAS-exon(80%A/T/G/C)

use strict; use warnings;
use READ_SAM;

print"usage: <polyAsite_finder_from_sam.pl> <wtts-aligned-sam-file>\n";
print"the out file is xxx.PAS.info\n";

my$FirstExonCutoff = 16;  #for the PAS on the first exon of a spliced mapping
my$SimplicityCut = 0.8;  #For the first exon

die "cmd err\n" unless @ARGV == 1;

open(IN,"$ARGV[0]");
open(OUT,">$ARGV[0].PAS.sam");
open(OUT2,">$ARGV[0].PAS.info");
print OUT2 "ref	dir	pos	count	MaxRead	cigar	sam-pos	sam-seq\n";
open(MULTIMAP,">$ARGV[0].PAS.2+");
open(VAIN,">$ARGV[0].NoPAS.sam");

my %read_collected;
my %pos_count;
my %pos2MaxRead;

while (<IN>){
	if(/^\@/){
        # skip headers
		print OUT ;
		print VAIN;
		next
	};

    # read sam
	chomp;
	my%sam = READ_SAM::General($_);
	my$direction;

	if($sam{flag} & 16 ){  # if reversely mapped
		$direction = "-";
		if(	READ_SAM::ContainPolyAsiteReverse(%sam) && 
			FirstExonSizeBigAndSeqComplex(
                $sam{cigar},$FirstExonCutoff,$sam{flag},$sam{seq})
		){
            # Call PAS
			$read_collected{$sam{qname}} ++;  # count PAS for this read
			print OUT "$_\n";  #print PAS.sam file 
			# print "contain PolyAsite:$sam{cigar}\t$sam{seq}\n$sam{qname}\n";
			
            if($read_collected{$sam{qname}} > 1){	#check for multi-location alignment
				print MULTIMAP "$sam{qname},$sam{flag},$sam{rname},$sam{pos},$sam{mapq},$sam{cigar},$sam{mrnm},$sam{mpos},$sam{isize},$sam{seq},$sam{qual},$sam{opt}";
            }
			else{
				my$position = READ_SAM::LastMapSite(%sam); # different from positive match
				my$location = $sam{rname}."\t".$direction."\t"."$position";
				$pos_count{$location}++;

				if(! exists $pos2MaxRead{$location}{length} || 
                    $pos2MaxRead{$location}{length} < $sam{cigar_length})
                { 
					$pos2MaxRead{$location}{length} = $sam{cigar_length}; 
					$pos2MaxRead{$location}{ReadID} = $sam{qname};
					$pos2MaxRead{$location}{cigar} = $sam{cigar};
					$pos2MaxRead{$location}{seq} = $sam{seq};
					$pos2MaxRead{$location}{mapq} = $sam{mapq};
					$pos2MaxRead{$location}{pos} = $sam{pos};
#					print "$pos2MaxRead{$location}{ReadID},$pos2MaxRead{$location}{length},$pos2MaxRead{$location}{cigar},$pos2MaxRead{$location}{seq}\n";
				}
			}
		}else{print VAIN "$_\n"}
	}
	else{
		$direction = "+";

		if(READ_SAM::ContainPolyAsite(%sam)  && 
            FirstExonSizeBigAndSeqComplex($sam{cigar},$FirstExonCutoff,$sam{flag},$sam{seq}))
        {
			print OUT "$_\n";#print sam file 
			$read_collected{$sam{qname}} ++;
	        #print "contain PolyAsite:$sam{cigar}\t$sam{seq}\n$sam{qname}\n";

			if($read_collected{$sam{qname}} > 1){
				print MULTIMAP "$sam{qname},$sam{flag},$sam{rname},$sam{pos},$sam{mapq},$sam{cigar},$sam{mrnm},$sam{mpos},$sam{isize},$sam{seq},$sam{qual},$sam{opt}\n";
            }
			else{
				my$position = $sam{pos};# the first map site
				my$location = $sam{rname}."\t".$direction."\t"."$position";
				$pos_count{$location}++;

				if(! exists $pos2MaxRead{$location}{length} || 
                    $pos2MaxRead{$location}{length} < $sam{cigar_length})
                { 
					$pos2MaxRead{$location}{length} = $sam{cigar_length}; 
					$pos2MaxRead{$location}{ReadID} = $sam{qname};
					$pos2MaxRead{$location}{cigar} = $sam{cigar};
					$pos2MaxRead{$location}{seq} = $sam{seq};
					$pos2MaxRead{$location}{mapq} = $sam{mapq};
					$pos2MaxRead{$location}{pos} = $sam{pos};
#					print "$pos2MaxRead{$location}{ReadID},$pos2MaxRead{$location}{length},$pos2MaxRead{$location}{cigar},$pos2MaxRead{$location}{seq}\n";
				}
			}
		}else{print VAIN "$_\n"}
	}
}

foreach (
    sort{
        (split /\t/,$a)[0] cmp (split /\t/,$b)[0] || 
        (split /\t/,$a)[1] cmp (split /\t/,$b)[1] || 
        (split /\t/,$a)[2] <=> (split /\t/,$b)[2] } 
        keys %pos_count)
{
	print OUT2 "$_	$pos_count{$_}	$pos2MaxRead{$_}{ReadID}	$pos2MaxRead{$_}{cigar}	$pos2MaxRead{$_}{pos}	$pos2MaxRead{$_}{mapq}	$pos2MaxRead{$_}{seq}\n"
}

print "finished//\n\n";

sub MD_QualifyPAS{
	# to detect the quality of the first exon
	my$md = shift;
	if($md =~ /MD:Z:(.*)/){$md = $1}
	else{print "$md is in wrong format\n"}
}

sub FirstExonSizeBigAndSeqComplex{
	#usage: my$test = FirstExonSizeBig($sam{cigar},$FirstExonCutoff,$sam{flag},$sam{seq});
	#function: 1.	First Exon Big enough to be specific
	#function: 2. The first exon can't be AA..AA/TT..TT/GG..GG/CC..CC/ such simple sequence
	my($cigar,$FirstExonCutoff,$flag,$SamSeq) = @_;
	my@exons = split/\d+N/,$cigar;
	#exon length and exon seq
	my%exon_length;my$ExonSeq;
	if($flag & 16){
		%exon_length = READ_SAM::CigarLength ($exons[-1]);
		$ExonSeq = substr $SamSeq,(length $SamSeq) - $exon_length{sam_seq_length};
	}else{
		%exon_length = READ_SAM::CigarLength ($exons[0]);
		$ExonSeq = substr $SamSeq,0,$exon_length{sam_seq_length};
	}

	my@ATGC= ("A","T","G","C","U");
	my$SimpleSeq = 0;
	foreach (@ATGC){
		my$bp_count = () = $ExonSeq =~ /($_)/gi;
		my$n_count = () = $ExonSeq =~ /(N)/gi;
		my$SimpleBp = $bp_count + $n_count;
		if($SimpleBp / length($ExonSeq) >= $SimplicityCut) {
			$SimpleSeq = 1 ; 
			last;
		}
	}

	if($exon_length{mapped_read_length} >= $FirstExonCutoff && $SimpleSeq == 0){return 1}
	else {return 0}
}
