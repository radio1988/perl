use strict; use warnings;

print "
usage: <pl> <xxx.fastq>

reads fastq and outputs fasta to std_out
";

my$count=0;
while(<>){
	$count++;
	if($count % 4 == 1) {$_ =~ s/^\@/>/;print}
	if($count % 4 == 2) {print}


}
