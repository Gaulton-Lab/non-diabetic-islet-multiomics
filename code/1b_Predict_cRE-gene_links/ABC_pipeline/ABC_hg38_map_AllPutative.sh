#!/usr/bin/env bash
### In this script we will map ALL ABC predictions (EnhancerPredictionsAllPutative.txt.gz) to hg38 by component.
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
gunzip -cd EnhancerPredictionsAllPutative.txt.gz | awk 'BEGIN{FS=OFS="\t"}{print $0,NR-1}'  > EnhancerPredictionsAllPutative.rownames.txt

#Remove the header and then perform liftover to hg38 for the CRE coordinations
tail -n+2 EnhancerPredictionsAllPutative.rownames.txt > EnhancerPredictionsAllPutative.rmHeader.txt
chain="/nfs/lab/ABC/references/hg19ToHg38.over.chain"
CrossMap.py region $chain EnhancerPredictionsAllPutative.rmHeader.txt EnhancerPredictionsAllPutative.crecoords.hg38.txt

#Make bed file of gene TSS coords from EnhancerPredictions and map to hg38
cut -f1,8 -d $'\t' EnhancerPredictionsAllPutative.rownames.txt > EnhancerPredictionsAllPutative.genecoords.prelim.txt
awk 'BEGIN{FS=OFS="\t"}{$3=$2+1; print $0; }' EnhancerPredictionsAllPutative.genecoords.prelim.txt > EnhancerPredictionsAllPutative.genecoords.prelim2.txt
awk 'BEGIN{FS=OFS="\t"}{print $0,NR-1}' EnhancerPredictionsAllPutative.genecoords.prelim2.txt > EnhancerPredictionsAllPutative.genecoords.prelim3.txt
tail -n+2 EnhancerPredictionsAllPutative.genecoords.prelim3.txt > EnhancerPredictionsAllPutative.genecoords.hg19.txt
CrossMap.py bed $chain EnhancerPredictionsAllPutative.genecoords.hg19.txt EnhancerPredictionsAllPutative.genecoords.hg38.txt

##### In R combine pieces so you can skip rows that didn't map #####
map_Rscript='/nfs/lab/ABC/code/ABC_hg38_map_AllPutative.R'
Rscript $map_Rscript $outdir


##### After finishing, check outputs and remove intermediates #####
#checking outputs
# head EnhancerPredictionsAllPutative.hg38.mapped.bedpe
# tail EnhancerPredictionsAllPutative.hg38.mapped.bedpe
wc -l EnhancerPredictionsAllPutative.hg38.mapped.bedpe
wc -l EnhancerPredictionsAllPutative.rmHeader.txt

rm EnhancerPredictionsAllPutative.rownames.txt
rm EnhancerPredictionsAllPutative.rmHeader.txt
rm EnhancerPredictionsAllPutative.crecoords*
rm EnhancerPredictionsAllPutative.genecoords*
