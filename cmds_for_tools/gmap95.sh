if [ "$#" -ne 2 ]; then
	echo "usage: <sh> <reads.fq> <threads>"
	exit 1
fi
READS=$1
THREADS=$2
gmap -d xt7 -f samse -n 0 --min-trimmed-coverage=0.9 --min-identity=0.95 -t $THREADS -B 4 --expand-offset=1 $READS > $READS.gmap95.sam 2> $READS.gmap95.sam.log
