#!/bin/sh

#Usage: sh call_featureCounts.sh [list of bam files to include] [saf file] [/path/to/output/file.txt]

bam_list=$( cat "$1" | tr '\r\n' ' ' )

echo "$bam_list"
featureCounts -p -T 10 -F SAF --donotsort -a "$2" -o "$3" $bam_list
