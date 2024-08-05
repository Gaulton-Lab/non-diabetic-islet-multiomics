#!/usr/bin/env bash
### In this script we will map the ABC predictions (EnhancerPredictions.txt) to hg38 by component.
### First we will separate the CRE and gene TSS coords and individually lift them over. Then we will
### combine them back together by only using rows that mapped in both processes.

### Read in dynamic inputs from flags using getopts
# Necessary input: directory containing the EnhancerPredictions.txt file (generally ABC_run/celltype/Prediction)
while getopts d: flag
do
    case "${flag}" in
        d) outdir=${OPTARG};;
    esac
done


##### Basic file manipulation with bash #####
#cd to the output dir with the EnhancerPredictions.txt to liftover
cd $outdir

#Add unique identifier to each line of EnhancerPredictions.txt (could be anything, easiest thing is just to add rownumbers starting at 0 as the col name)
awk 'BEGIN{FS=OFS="\t"}{print $0,NR-1}' EnhancerPredictions.txt  > EnhancerPredictions.rownames.txt

#Remove the header and then perform liftover to hg38 for the CRE coordinations
tail -n+2 EnhancerPredictions.rownames.txt > EnhancerPredictions.rmHeader.txt
chain="/nfs/lab/ABC/references/hg19ToHg38.over.chain"
CrossMap.py region $chain EnhancerPredictions.rmHeader.txt EnhancerPredictions.crecoords.hg38.txt

#Make bed file of gene TSS coords from EnhancerPredictions and map to hg38
cut -f1,6 -d $'\t' EnhancerPredictions.rownames.txt > EnhancerPredictions.genecoords.prelim.txt
awk 'BEGIN{FS=OFS="\t"}{$3=$2+1; print $0; }' EnhancerPredictions.genecoords.prelim.txt > EnhancerPredictions.genecoords.prelim2.txt
awk 'BEGIN{FS=OFS="\t"}{print $0,NR-1}' EnhancerPredictions.genecoords.prelim2.txt > EnhancerPredictions.genecoords.prelim3.txt
tail -n+2 EnhancerPredictions.genecoords.prelim3.txt > EnhancerPredictions.genecoords.hg19.txt
CrossMap.py bed $chain EnhancerPredictions.genecoords.hg19.txt EnhancerPredictions.genecoords.hg38.txt


##### In R combine pieces so you can skip rows that didn't map #####
map_Rscript='/nfs/lab/ABC/code/ABC_hg38_map.R'
Rscript $map_Rscript $outdir


##### After finishing, check outputs and remove intermediates #####
#checking outputs
# head EnhancerPredictions.hg38.mapped.bedpe
# tail EnhancerPredictions.hg38.mapped.bedpe
wc -l EnhancerPredictions.hg38.mapped.bedpe
wc -l EnhancerPredictions.rmHeader.txt

#remove intermediates
rm EnhancerPredictions.rownames.txt
rm EnhancerPredictions.rmHeader.txt
rm EnhancerPredictions.crecoords*
rm EnhancerPredictions.genecoords*
