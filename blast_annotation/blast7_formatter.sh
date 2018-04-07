#!/bin/bash
echo "Function:	transform blast result from blast11(asn1) to blast7(tab.txt)"
echo "Usage:	<blast&_formatter.sh> <*.blast11>"
blast11_file=$1
echo "working in the background"
blast_formatter -archive $blast11_file -out $blast11_file.blast7 -outfmt '7 qseqid qlen qcovs pident sseqid slen length nident positive staxids sscinames stitle evalue bitscore sstrand'

#nohup blast_formatter -archive $blast11_file -out $blast11_file.Blast7 -outfmt '7 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp' &> $blast11_file.Blast7.log &
#echo 'nohup blast_formatter -archive $blast11_file -out $blast11_file.Blast7 -outfmt '7 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp' &> $blast11_file.Blast7.log &' > $blast11_file.Blast7.cmd
#44 parameters in total

