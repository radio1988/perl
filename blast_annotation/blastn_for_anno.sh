echo 'parameters optimised for blastn for remote species'
echo 'usage : <sh> <query.fa> <threads>'
in=$1
threads=$2
echo'
blastn -task blastn -query $in -db refseq_rna -gilist ~/db/blast/gil/7711.chordata.ref_nucl.gil -evalue 0.001 -outfmt 11 -num_threads $threads -max_target_seqs 100 -out $in.7711_refrna.asn
'
nohup blastn -task blastn -query $in -db refseq_rna -gilist ~/db/blast/gil/7711.chordata.ref_nucl.gil -evalue 0.000001 -outfmt 11 -num_threads $threads -max_target_seqs 10 -out $in.7711_refrna.asn &> $in.7711_refrna.log

