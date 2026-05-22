#!/bin/bash

r1=$1
r2=$2
out_r1=$3
out_r2=$4
build=$5
threads=${6:-8}

if [ $# -lt 5 ];then
	echo "Not enough arguments. Usage: $0 <R1> <R2> <out_R1> <out_R2> <build> <threads>"
	exit 1
fi

bbmap.sh build=$build in=$r1 in2=$r2 outm=${r1/.fastq.gz/}.tmp-bbmap.fastq minid=0.5 threads=$threads
reformat.sh in=${r1/.fastq.gz/}.tmp-bbmap.fastq out1=${out_r1/.gz/} out2=${out_r2/.gz/}
rm ${r1/.fastq.gz/}.tmp-bbmap.fastq

if [[ "$(echo ${out_r1##*.})" == "gz" ]] && [[ "$(echo ${out_r2##*.})" == "gz" ]];then 
    gzip ${out_r1/.gz/} ${out_r2/.gz/}
fi
