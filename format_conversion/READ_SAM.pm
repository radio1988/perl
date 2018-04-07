package READ_SAM;
use strict;use warnings;
sub General{
	#only parse one line
	my $line = shift;
	$line =~ s/[\r\n]//g;
	if($line =~ /^\@/) {return 0};
	my %sam;
	($sam{qname},$sam{flag},$sam{rname},$sam{pos},$sam{mapq},$sam{cigar},$sam{mrnm},$sam{mpos},$sam{isize},$sam{seq},$sam{qual},$sam{opt}) = split /\t/,$line, 12;
	if($sam{flag} == 4) {return %sam}

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
	$sam{cigar_length} = $scount + $hcount + $mcount + $icount + $Xcount + $Equal_count;#length of the read to find the longest read for cuffcompare annotation
	$sam{last_mapp_site} = ( $sam{pos} - 1 ) + $mcount + $dcount + $ncount + $Xcount + $Equal_count;

	($sam{NM}) = $sam{opt} =~ /NM:i:(\d+)/;
	$sam{mapped_read_length} = $mcount + $icount + $Xcount + $Equal_count; #length of mapped proportion of the read
	$sam{clipped_identity} = 1 - $sam{NM}/$sam{mapped_read_length};
	$sam{coverage} = $sam{mapped_read_length}/$sam{cigar_length};

	#sam mini_anchoring length
	my@exon_cigar = split /N/,$sam{cigar};
	my$anchor_exon1 = $exon_cigar[0];
	my$anchor_exon2 = $exon_cigar[-1];
	#$sam{min_exon_length} = 

	#GC content, consider only ATGC, U not allowed, N ignored
	my$GC_count = () = $sam{seq} =~ /([GC])/g;
	my$AT_count = () = $sam{seq} =~ /([AT])/g;
	$sam{GC_content} = $GC_count/($GC_count+$AT_count);

	if($sam{flag} & 16){$sam{dir} = "-"}else{$sam{dir} = "+"}
	
	return (%sam);
}

sub Print {
	my(%sam) = @_;
	return ("$sam{qname}	$sam{flag}	$sam{rname}	$sam{pos}	$sam{mapq}	$sam{cigar}	$sam{mrnm}	$sam{mpos}	$sam{isize}	$sam{seq}	$sam{qual}	$sam{opt}")

}

sub CigarLength{
	my%cigar;
	my$cigar = shift;
	my@cigar_S = $cigar =~ /(\d+)S/g;
        my@cigar_H = $cigar =~ /(\d+)H/g;
        my@cigar_M = $cigar =~ /(\d+)M/g;
        my@cigar_I = $cigar =~ /(\d+)I/g;
        my@cigar_D = $cigar =~ /(\d+)D/g;
        my@cigar_N = $cigar =~ /(\d+)N/g;
        my@cigar_X = $cigar =~ /(\d+)X/g;
        my@cigar_euqal = $cigar =~ /(\d+)=/g;
        my($mcount,$icount,$dcount,$scount,$hcount,$ncount,$Xcount,$Equal_count) = (0,0,0,0,0,0,0,0);
        foreach(@cigar_S){$scount += $_}
        foreach(@cigar_H){$hcount += $_}
        foreach(@cigar_M){$mcount += $_}
        foreach(@cigar_I){$icount += $_}
        foreach(@cigar_D){$dcount += $_}
        foreach(@cigar_N){$ncount += $_}
        foreach(@cigar_X){$Xcount += $_}
        foreach(@cigar_euqal){$Equal_count += $_}
        $cigar{cigar_length} = $scount + $hcount + $mcount + $icount + $Xcount + $Equal_count;#length of the read to find the longest read for cuffcompare annotation
        $cigar{sam_seq_length} = $scount + $mcount + $icount + $Xcount + $Equal_count;#length of the sequence shown in sam file
        $cigar{mapped_read_length} = $mcount + $icount + $Xcount + $Equal_count; #length of mapped proportion of the read
	$cigar{intron_length} = $ncount;
	return %cigar;
}

sub ContainPolyAsite{
	my(%sam) = @_;#a hash of sam single line
#	print "test:$sam{qname};\n";	
	my($soft_clip,$hard_clip) ;
	my$max_hang = 5;
#	for CTTT header, clip CTTT in counting
	if ($sam{cigar} =~ /^(\d+)S/){$soft_clip = $1}	
	if ($sam{cigar} =~ /^(\d+)H/){$hard_clip = $1}	

#	if the mapping point is within 5bp of the read, use it, other wise discard it
	if (	$sam{cigar} =~ /^\d+M/) {return 1}
	elsif(	$sam{cigar} =~ /^\d+S/ && $soft_clip <= $max_hang){return 1}
	elsif( $sam{cigar} =~ /^\d+H/ &&  $hard_clip <= $max_hang ){return 1}

	return 0;
}

sub ContainPolyAsiteReverse{
	my(%sam) = @_;#a hash of sam single line
#	print "test:$sam{qname};\n";	
	my($soft_clip,$hard_clip) ;
	my$max_hang = 5;
#	for CTTT header, clip CTTT in counting
	if ($sam{cigar} =~ /(\d+)S$/){$soft_clip = $1}	
	if ($sam{cigar} =~ /(\d+)H$/){$hard_clip = $1}	

#	if the mapping point is within 5bp of the read, use it, other wise discard it
	if (	$sam{cigar} =~ /\d+M$/) {return 1}
	elsif(	$sam{cigar} =~ /\d+S$/ && $soft_clip <= $max_hang){return 1}
	elsif( $sam{cigar} =~ /\d+H$/ &&  $hard_clip <= $max_hang ){return 1}

	return 0;
}

sub LastMapSite{
	my(%sam) = @_;
	my@cigar_s = $sam{cigar} =~ /(\d+)S/g;	
	my@cigar_H = $sam{cigar} =~ /(\d+)H/g;	
	my@cigar_M = $sam{cigar} =~ /(\d+)M/g;	
	my@cigar_I = $sam{cigar} =~ /(\d+)I/g;	
	my@cigar_D = $sam{cigar} =~ /(\d+)D/g;	
	my@cigar_N = $sam{cigar} =~ /(\d+)N/g;	
	my@cigar_X = $sam{cigar} =~ /(\d+)X/g;	
	my@cigar_euqal = $sam{cigar} =~ /(\d+)=/g;	
	#foreach(@cigar_M) {print "$_\t"}
	#print "\n";
	my($mcount,$icount,$dcount,$scount,$hcount,$ncount,$Xcount,$Equal_count) = (0,0,0,0,0,0,0,0);
	foreach(@cigar_M){$mcount += $_}
	foreach(@cigar_I){$icount += $_}
	foreach(@cigar_D){$dcount += $_}
	foreach(@cigar_N){$ncount += $_}
	foreach(@cigar_X){$Xcount += $_}
	foreach(@cigar_euqal){$Equal_count += $_}
	my$polyAsite = ( $sam{pos} - 1 ) + $mcount + $dcount + $ncount + $Xcount + $Equal_count;
	#print"$sam{pos},$sam{cigar},$polyAsite\n\n";
	return ($polyAsite);
}

1;

