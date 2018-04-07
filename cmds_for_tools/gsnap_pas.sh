echo "usage: <sh> <fq> <threads>"
if [ "$#" -ne 2 ]; then
        echo "command error"
        exit 1
fi



READ=$1
TREAD=$2

gsnap -d xt7 -s 7.2 -N 1 -m 0.05 -t $TREAD -A sam -n 1 -B 5 --failed-input=$READ.unmapped $READ > $READ.gsnap.sam
