#!/usr/bin/perl
#Rui Li, Animal Sciences Department, Washington State University;
#PhD student of Dr. Zhihua Jiang
#Email: liruiradiant+perl@gmail.com

use strict;
use warnings;
use READ_SYMBOL;

##input:
#*infos.JGI_symbol:	(scaffold_1      +       35441   2       EN03R:01845:13911       5S107M6S        35441   40      GTTTTTGTATAACAATTGTGTATTTAATATCACATAAAACACATTCAGTACTTAAAACAAGCGACGTCGATCAAATGTATGTACAAAGCTTTACATGTGTTTGAAGTACAAAAGGCCC  TTAAGTACTGAATGTGTTTTATGTGATATTAAATACACAATTGTTATACAgatctgccttgtgtgggttatttaacttgctaacgtttgtctttattgcc       Xetro.A00001    c       EN03R:01845:13911       EN03R:01845:13911       Xetro.A00001    XB-GENE-488014  aqp3    aquaporin 3 (Gill blood group))

#the input file don't have to be sorted
#there can be a header of file(ref	dir)

##output: *clust/clust2/clust.err(for multiple gene_acc)

##function 
#1.single-linkage clustering
#2.output raw clusters into *clust2
#3.calculate class_code
#4.calculate the expression of each cluster
#5.calculate the average location
#6.output summarised clusters into *clust

#this is clustering version2,input is *JGI_symbol
#2014-11-5 bug of chromosome change fixed
#2014-11-17 add annotation info to *clust by reading *infos.JGI_symbol file
#2014-12-18 add max-read and 100seq_of_average_site into *clust
#2015-1-5 fix ID bug in *clust (wrongly collected mapped ID in clolumn 15, should be column 11)
#2015-1-5 fix class_code 'u' bug ( with ID but code is 'u')


#note
#rname: "scaffold_x"
#direction: "+/-"
#location: "10000"
#count_site: "20" (the count for each 1bp precise site)
my $usage = "usage: <popl> <*infos.JGI_symbol> <cutoff>\n";
die "command err: $usage\n" unless @ARGV == 2;
print"> polyAsite_clustering.pl $ARGV[0] $ARGV[1]\n";
my$CutOff = $ARGV[1];
open(IN,$ARGV[0]) or die "reading info file err\n";
open(OUT,">$ARGV[0].$ARGV[1]clust") or die "err output into file\n";
open(OUT2,">$ARGV[0].$ARGV[1]clust2") or die "err output into file2\n";
open(ERR,">$ARGV[0].$ARGV[1]clust.err") or die "err output into err file\n";

my($rname_old,$direction_old,$pos_old,$count_old);
my$line_n = 0;
my@cluster_t;

my@in = <IN>;
my@clust2;

my@clust;#one summary line for each cluster

###clustering
my$i=0;#for indexing the clusters
foreach( 
		#sort by scaffold_id, direction and location
		sort{
		(split /\t/,$a)[0] cmp (split /\t/,$b)[0] ||
		(split /\t/,$a)[1] cmp (split /\t/,$b)[1] ||
		(split /\t/,$a)[2] <=> (split /\t/,$b)[2] 
		}  @in
       )
	{
		#header
		if(/^ref\tdir/){next}
		#prepare
		$line_n++;
		chomp;
		#collect location info:	#scaffold_1      +       35441   2 
		if($line_n == 1) { ($rname_old,$direction_old,$pos_old,$count_old) = split /\t/,$_ }

		my($rname,$direction,$pos,$count)= split /\t/,$_;
		if(
			WithinRange($rname,$direction,$pos,$count,$rname_old,$direction_old,$pos_old,$count_old)
		)
		{
			push @cluster_t, $_;
		}
		else{	
			#not within range, start a new cluster
			$i ++;	#cluster number ++
			my $j=0;#expression level
			foreach(@cluster_t){
				$j++;
				$clust2[$i][$j] = $_;#store result in a 2D-array
			}

			my($mean_rname,$mean_dir,$mean_pos,$mini_pos,$max_pos,$sum_count) = ClusterOut(@cluster_t);
			my$out1 = "$mean_rname\t$mean_dir\t$mean_pos\t$mini_pos-$max_pos\t$sum_count";
			$clust[$i] = $out1; 

			#reset for next round
			@cluster_t = ();
			push @cluster_t, $_;
		}
		($rname_old,$direction_old,$pos_old,$count_old) = ($rname,$direction,$pos,$count);
	}

#fix the last cluster
$i ++;	
my$j=0;
foreach(@cluster_t){
	$j++;
	$clust2[$i][$j] = $_;
}
my($mean_rname,$mean_dir,$mean_pos,$mini_pos,$max_pos,$sum_count) = ClusterOut(@cluster_t);
my$out1 = "$mean_rname\t$mean_dir\t$mean_pos\t$mini_pos-$max_pos\t$sum_count";
$clust[$i] = $out1; 

###OUT2 *clust2
print OUT2 "cluster_index\tscaffold\tdirection\tposition\trange\tcount\tclass_code\tgene_acc\tgene_page\tgene_symbol\tgene_desc\n";
for my$i ( 1 .. $#clust2){
	print OUT2 "#cluster$i:$clust[$i]\n";
	for my$j ( 1 .. $#{$clust2[$i]}) {
		print OUT2 "$clust2[$i][$j];\n";
	}
}

###OUT1 *clust
for my$i ( 1 .. $#clust2){
	print OUT "#cluster$i:\t";
	my(%gene_acc,%ac_anno,%class_code,@read);
	my$class_code="";

	#prepare
	for my$j ( 1 .. $#{$clust2[$i]}) {
		my%pas = READ_SYMBOL::SYMBOL($clust2[$i][$j]);
		if($pas{gene_acc} ne "-"){ 
			$gene_acc{$pas{gene_acc}} ++;
			$ac_anno{$pas{gene_acc}} = "-";
		}
		if($pas{class_code} ne "-") {$class_code{$pas{class_code}} ++;}
		#for the max-read infomation
		$read[$j]{ID} = $pas{read_id};
		$read[$j]{seq} = $pas{read_seq};
		$read[$j]{cigar} = $pas{cigar};
		$read[$j]{mq} = $pas{mq};
		
	}

	#prepare for class-code
	my($max_class_code,$max_count)= ("",0);#the class code with largest count
	($class_code,$max_class_code,$max_count) = CLASS_CODE(%class_code);

	#prepare for max-read
	my($max_read_id,$max_read_cigar,$max_read_seq,$max_read_mq,$max_read_length)= ("","","","",0);
	foreach my$j(1..$#read){
		if (length $read[$j]{seq} > $max_read_length){
			$max_read_id = $read[$j]{ID};
			$max_read_seq = $read[$j]{seq};
			$max_read_mq = $read[$j]{mq};
			$max_read_cigar = $read[$j]{cigar};
			$max_read_length = length $read[$j]{seq};
		}
	}

	#out
	if( scalar keys %gene_acc > 1){
		print ERR "#cluster$i:", scalar keys %gene_acc,"\n";#number of accessions
		foreach (keys %gene_acc){print ERR "$_	$gene_acc{$_}	$ac_anno{$_};\n"}#each accession
		print ERR "\n";

		print OUT "$clust[$i]\tERR","\t-"x2,"\t$max_read_id\t$max_read_mq\t$max_read_cigar\t$max_read_length\t$max_read_seq\n";#mark for manual annotation, because cuffcompare returned multiple gene_accession

	}elsif( scalar keys %gene_acc == 0){
		print OUT "$clust[$i]\tu","\t-"x2,"\t$max_read_id\t$max_read_mq\t$max_read_cigar\t$max_read_length\t$max_read_seq\n";#may cause bug for non-standard cuffcompare input?

	}else{
		my$gene_acc = (keys %gene_acc)[0];
		print OUT "$clust[$i]\t$max_class_code:$class_code\t$gene_acc\t$ac_anno{$gene_acc}\t$max_read_id\t$max_read_mq\t$max_read_cigar\t$max_read_length\t$max_read_seq\n";

	}
}

#subs
sub WithinRange{#if within range, return 1
	my($rname,$direction,$pos,$count,$rname_old,$direction_old,$pos_old,$count_old) = @_;
	if( 		$rname ne $rname_old ||
			$direction_old ne $direction ||
			abs($pos - $pos_old) > $CutOff
	  ){return 0}

	return 1;

}

sub CLASS_CODE{
	my(%class_code) = @_;	
	my($class_code,$max_class_code,$max_count) = ("","",0);
	my@class_sort;
	foreach(sort { $class_code{$b} <=> $class_code{$a} } keys %class_code){
		$class_code .= "$class_code{$_}$_/";
		push @class_sort, "$_\t$class_code{$_}" ;
	}
	
	if(scalar @class_sort > 1){
		if($class_sort[0] =~ /u\t/){
			($max_class_code,$max_count) = split "\t",$class_sort[1];
		}else{
			($max_class_code,$max_count) = split "\t",$class_sort[0];
		} 
	}else{
		($max_class_code,$max_count) = split "\t",$class_sort[0];
	}

	return ($class_code,$max_class_code,$max_count);
}

sub ClusterOut{
	my@cluster_t = @_;
	my$sum_count;
	my$sum_pos_count;
	my($rname,$direction,$pos,$count); 
	foreach(@cluster_t){
		chomp;
		($rname,$direction,$pos,$count) = split /\t/,$_;
		$sum_count += $count;
		$sum_pos_count += ($count * $pos);
	}
	my$mini_pos = (split/\t/,$cluster_t[0])[2];
	my$max_pos = (split/\t/,$cluster_t[-1])[2];
	my$mean_pos = int($sum_pos_count/$sum_count+0.5);
	return($rname,$direction,$mean_pos,$mini_pos,$max_pos,$sum_count);

}
