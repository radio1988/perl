package READ_GFF;
use strict;use warnings;
sub CUFFLINKS{
	my$line = shift;
	chomp $line;
	my@ele = split /\t/,$line;
	my %gff ;
	(
		$gff{seqid},
		$gff{source},
		$gff{type},
		$gff{start},
		$gff{end},
		$gff{score},
		$gff{strand},
		$gff{phase},
		$gff{attributes}
	)=@ele;
	unless ( defined $gff{attributes} ) {die "format err in $line;\n"};	

	if($gff{attributes} =~ /gene_id "([^"]*)"/){ $gff{gene_id} = $1;}
	elsif( $gff{attributes} =~ /Name=([^;]*);/){ $gff{gene_id} = $1;}

	if($gff{attributes} =~ /transcript_id "([^"]*)"/){ $gff{transcript_id} = $1}
	if($gff{attributes} =~ /gene_name "([^"]*)"/){ $gff{gene_name} = $1}
	if($gff{attributes} =~ /exon_number "([^"]*)"/){ $gff{exon_number} = $1	}
	if($gff{attributes} =~ /class_code "([^"]*)"/){ $gff{class_code} = $1}
	if($gff{attributes} =~ /nearest_ref "([^"]*)"/){ $gff{nearest_ref} = $1}

	return %gff;


}
sub TransDecoder{
	my$line = shift;
	chomp $line;
	$line = $line.";" ;
	$line =~ s/ //g;
	my@ele = split /\t/,$line;
	my %gff ;
	(
		$gff{seqid},
		$gff{source},
		$gff{type},
		$gff{start},
		$gff{end},
		$gff{score},
		$gff{strand},
		$gff{phase},
		$gff{attributes}
	)=@ele;
	unless ( defined $gff{attributes} ) {die "format err in $line;\n"};	

	if( $gff{attributes} =~ /ID=([^;]*);/){ $gff{ID} = $1;}
	if( $gff{attributes} =~ /Parent=([^;]*);/){ $gff{parent} = $1;}
	if( $gff{attributes} =~ /geneID=([^;]*);/){ $gff{geneID} = $1;}
	if( $gff{attributes} =~ /gene_name=([^;]*);/){ $gff{gene_name} = $1;}

	return %gff;
}
sub XENBASE{
	my$line = shift;
	chomp $line;
	my@ele = split /\t/,$line;
	my %gff ;
	(
		$gff{seqid},
		$gff{source},
		$gff{type},
		$gff{start},
		$gff{end},
		$gff{score},
		$gff{strand},
		$gff{phase},
		$gff{attributes}
	)=@ele;
	
	if($gff{attributes} =~ /Name=([^;]*);/){ $gff{symbol} = $1;}
	if($gff{attributes} =~ /Name=([^;]*);/){ $gff{name} = $1;}
	if($gff{attributes} =~ /Name=([^;]*);/){ 
		$gff{JGI_ID} = $1;
		($gff{JGI_ID}) = $gff{JGI_ID} =~ /(Xetro\.[^\.]*)\./;
	}
	if($gff{attributes} =~ /ID=([^;]*);/){ $gff{ID} = $1;}
	if($gff{attributes} =~ /ID=([^;]*);/){ $gff{PAC_ID} = $1;}
	if($gff{attributes} =~ /pacid=([^;]*);/){ $gff{pacid} = $1;}
	if($gff{attributes} =~ /Parent=([^;]*);/){ $gff{parent} = $1;}

	return %gff;


}
sub MITO{
	die"not finished 20160209";
	my$line = shift;
	chomp $line;
	my@ele = split /\t/,$line;
	my %gff ;
	(
		$gff{seqid},
		$gff{source},
		$gff{type},
		$gff{start},
		$gff{end},
		$gff{score},
		$gff{strand},
		$gff{phase},
		$gff{attributes}
	)=@ele;
	
	if($gff{attributes} =~ /Name=([^;]*);/){ $gff{symbol} = $1;}
	if($gff{attributes} =~ /Name=([^;]*);/){ $gff{name} = $1;}
	if($gff{attributes} =~ /Name=([^;]*);/){ 
		$gff{JGI_ID} = $1;
		($gff{JGI_ID}) = $gff{JGI_ID} =~ /(Xetro\.[^\.]*)\./;
	}
	if($gff{attributes} =~ /ID=([^;]*);/){ $gff{ID} = $1;}
	if($gff{attributes} =~ /ID=([^;]*);/){ $gff{PAC_ID} = $1;}
	if($gff{attributes} =~ /pacid=([^;]*);/){ $gff{pacid} = $1;}
	if($gff{attributes} =~ /Parent=([^;]*);/){ $gff{parent} = $1;}

	return %gff;


}
sub GFF121{ #GCF_000471725.1_UMD_CASPUR_WB_2.0_genomic.gff
	my$line = shift;
	chomp $line;
	my@ele = split /\t/,$line;
	my %gff ;
	(
		$gff{seqid},
		$gff{source},
		$gff{type},
		$gff{start},
		$gff{end},
		$gff{score},
		$gff{strand},
		$gff{phase},
		$gff{attributes}
	)=@ele;

	if($gff{attributes} =~ /ID=([^;]*);/){ $gff{ID} = $1;}
	if($gff{attributes} =~ /GeneID:([^;]*);/){ $gff{Gene_ID} = $1;}
	if($gff{attributes} =~ /Genbank:([^;]*);/){ $gff{Genbank} = $1;}
	if($gff{attributes} =~ /Name=([^;]*);/){ $gff{name} = $1;}
	if($gff{attributes} =~ /gbkey=([^;]*);/){ $gff{gbkey} = $1;}
	if($gff{attributes} =~ /mol_type=([^;]*);/){ $gff{mol_type} = $1;}
	if($gff{attributes} =~ /Parent=([^;]*);/){ $gff{parent} = $1;}
	if($gff{attributes} =~ /gene=([^;]*);/){ $gff{gene} = $1;}
	if($gff{attributes} =~ /transcript_id=([^;]*)$/){ $gff{transcript_id} = $1;}
	if($gff{attributes} =~ /gene_biotype=([^;]*)$/){ $gff{gene_biotype} = $1;}
	
	
	if($gff{attributes} =~ /pacid=([^;]*);/){ $gff{pacid} = $1;}

	return %gff;


}

sub HASH1{
	#works only for  Xentr7_2_Stable.gff3 
	#input @gff (<GFF>)
	my@gff = @_;@_ = ();#the whole gff file

	my%JGI2mRNA;
	my%mRNA2CDS;

	foreach(@gff){
		$_ =~ s/[\r\n]//g;
		if(/^#/ || length($_) < 2){next}

		my%line = XENBASE($_);
		#JGI2mRNA
		if ($line{type} eq "mRNA"){
#			print "$_\n";
#			print "JGI: $line{JGI_ID};\nmRNA:$line{PAC_ID};\n\n";
			push @{ $JGI2mRNA{$line{JGI_ID}} }, $line{PAC_ID}; 

#			if(scalar @{ $JGI2mRNA{$line{JGI_ID}} } > 3){
#				print "scalar",scalar @{ $JGI2mRNA{$line{JGI_ID}} },"\n\n";
#				print "array:\n @{ $JGI2mRNA{$line{JGI_ID}} }\\n\n";	
#	#			sleep 1;
#			}
		}
		#mRNA2CDS
		if ($line{type} eq "CDS") {
#			print "CDS:$_\n";
			push @{ $mRNA2CDS{$line{parent}} }, $_; 
#			if(scalar @{ $mRNA2CDS{$line{parent}} } > 2){
#				foreach my$i (0 .. $#{ $mRNA2CDS{$line{parent}} }){
#					print"$i:   $mRNA2CDS{$line{parent}}[$i]\n";
#				}
#				print "\n\n";
#				sleep 2;
#			}
		}
		
	
	}

	return(\%JGI2mRNA,\%mRNA2CDS);
}


1;
