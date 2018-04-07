#!/bin/bash
echo 'makeblastdb -in $fastafile -input_type $intype -dbtype $dbtype'
echo "input filename"
read fastafile
echo "in_typei: fasta, blastdb,asn1_txt, asn1_bin"
read intype
echo "dbtype: prot, nucl"
read dbtype
makeblastdb -in	$fastafile -input_type $intype -dbtype $dbtype
