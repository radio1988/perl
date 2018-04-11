use strict; use warnings;

print "
usage: <pl> <xxx.fastq> <outname.fasta>
reads fastq and outputs fasta to std_out
";

die "cmd err\n" unless @ARGV==2;

open(IN, "$ARGV[0]") or die 'err reading\n';
open(OUT, ">$ARGV[1]") or die 'err writing\n';

my$count=0;
while(<IN>){
	$count++;
	if($count % 4 == 1) {$_ =~ s/^\@/>/;print OUT}
	if($count % 4 == 2) {print OUT}
}

print "conversion finished\n";

close IN;
close OUT;