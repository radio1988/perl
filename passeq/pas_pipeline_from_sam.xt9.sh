#!/bin/bash
echo "
>> <pas_pipeline.sh> <aligner.SAM> <PREFIX> 
<PREFIX> should be dynamic, containing sample information in a shell loop\n"
#201411: filter MQ, filter identity95, filter NNN, filter IP(A rich region PAS), auto get bam file
#20141218: added 100bp seq for *clust file
#20150126: add Mapped-length filter in sam_pas_filter.pl, no need to filter read length in trimming step(ttt50bp->ttt1bp)
#20150131: try not running two pipe in same folder for same fastq file at the same time
#20160212: start from filtered sam
if [ "$#" -ne 2	]; then
	echo "command error in pas_pipeline.sh"
	exit 1
fi

SAM=$1
#MQCUT=$2
PREFIX=$2

GENOME="/home/rui/db/genome9.0/add_mito/Xtropicalis.v9.mito.fa"
REFGFF="/home/rui/db/genome9.0/xt9_anno/Xtropicalisv9.0.Named.primaryTrs.gff3"
#GENE_LIST="$HOME/db/genome7.1/trans/stable/gene.sort.map"
#IDENTITYCUT=0.95
#COVERAGE=0.9
#MAPPED_LENGTH=16
echo $GENOME
echo $REFGFF
echo "CMD:pas_pipeline.sh $SAM $MQCUT $PREFIX 
PREFIX:$PREFIX
"
sleep 3

#filter sam
#sam_pas_quality_filter.pl $SAM $MQCUT $IDENTITYCUT $COVERAGE $MAPPED_LENGTH
#SAM=$SAM.hq$MQCUT.$IDENTITYCUT.$COVERAGE.$MAPPED_LENGTH

#find PAS from sam
polyAsite_finder_from_sam.pl $SAM
samtools_sam2bam.pl $SAM.PAS.sam &
polyAsite_100bpseq.pl $SAM.PAS.info  $GENOME
polyAsite_NNN_filter.pl $SAM.PAS.infos

#use cuffcompare for annotation
polyAsiteInfoCigar2gtf.pl $SAM.PAS.infos.N
gff_reverse_direction.pl $SAM.PAS.infos.N.gtf > $SAM.PAS.infos.N.rev.gtf
cuffcompare -G -r $REFGFF $SAM.PAS.infos.N.rev.gtf -o $PREFIX
cut -f 1,3,4,5  $PREFIX.$SAM.PAS.infos.N.rev.gtf.tmap >  $PREFIX.$SAM.PAS.infos.N.rev.gtf.tmap.txt

#add Gene symbol based on XenbaseID
Pmatch-mysql.pl $SAM.PAS.infos.N 5 $PREFIX.$SAM.PAS.infos.N.rev.gtf.tmap.txt 3
mv $SAM.PAS.infos.N.Pmatched $PREFIX.anno

#filter internal priming
polyAsite_internal_priming_filter.pl $PREFIX.anno #out put .ip .pass

#clustering
polyAsite_clustering.pl $PREFIX.anno.pass 5
polyAsite_clustering.pl $PREFIX.anno.ip 5
polyAsite_clustering.pl $PREFIX.anno 5
polyAsite_100bpseq.pl $PREFIX.anno.pass.5clust  $GENOME &
polyAsite_100bpseq.pl $PREFIX.anno.ip.5clust  $GENOME & 
polyAsite_100bpseq.pl $PREFIX.anno.5clust  $GENOME

polyAsite_gene_level_expression.pl $PREFIX.anno.pass.5clusts 8 6 7
polyAsite_gene_level_expression.pl $PREFIX.anno.ip.5clusts 8 6 7
polyAsite_gene_level_expression.pl $PREFIX.anno.5clusts 8 6 7

polyAsite_cluster_stat.pl $PREFIX.anno.pass.5clusts > $PREFIX.anno.pass.5clusts.stat
polyAsite_cluster_stat.pl $PREFIX.anno.ip.5clusts > $PREFIX.anno.ip.5clusts.stat
polyAsite_cluster_stat.pl $PREFIX.anno.5clusts > $PREFIX.anno.5clusts.stat

#get GE2 table
#Pmatch-mysql.pl $GENE_LIST 1 $PREFIX.anno.pass.5clusts.GE 1 && mv $GENE_LIST.Pmatched $PREFIX.anno.pass.5clusts.GE2
#Pmatch-mysql.pl $GENE_LIST 1  $PREFIX.anno.ip.5clusts.GE 1 && mv $GENE_LIST.Pmatched $PREFIX.anno.ip.5clusts.GE2
#Pmatch-mysql.pl $GENE_LIST 1  $PREFIX.anno.5clusts.GE 1 && mv $GENE_LIST.Pmatched $PREFIX.anno.5clusts.GE2


#clean unimportant data
rm -rf ${PREFIX}_extra_data
mkdir ${PREFIX}_extra_data
mv $PREFIX*10clust* ${PREFIX}_extra_data
mv $PREFIX*30clust* ${PREFIX}_extra_data
mv $PREFIX*clust ${PREFIX}_extra_data
mv $PREFIX*tracking ${PREFIX}_extra_data
mv $PREFIX*clust2 ${PREFIX}_extra_data
mv $PREFIX*tmap ${PREFIX}_extra_data
mv $PREFIX*refmap ${PREFIX}_extra_data
mv $PREFIX*loci ${PREFIX}_extra_data
mv $PREFIX*err ${PREFIX}_extra_data

mv $PREFIX*gtf ${PREFIX}_extra_data
#delete some temp data
#rm -f $PREFIX*JGI_ID $PREFIX*info






