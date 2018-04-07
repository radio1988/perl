if [ "$#" -ne 2 ]; then
        echo "command error: usage: <gsnap_xt9.sh> <*fq> <threads>"
        exit 1
fi

READ=$1
TREAD=$2
echo "gsnap_xt9.sh $READ $TREAD"

gsnap -d xt9.mito.k15 -s xt9.both.iit -N 1 --gunzip -t $TREAD -A sam -B 5 --failed-input=$READ.unmapped $READ > $READ.gsnap.sam 2> $READ.gsnap.sam.log
samtools_sam2bam.pl $READ.gsnap.sam

#-m 0.05
