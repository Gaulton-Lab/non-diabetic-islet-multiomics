#!/usr/bin/env bash
### In this script we will call peaks in hg19 for use in running ABC. This assumes we have the celltype specific tagAlign
### files from calling peaks on the same group of cell types that we plan to run ABC on!
### This assumes you have a conda env with macs2 installed. You can copy this one if necessary:
### conda create --name call_peaks --clone /home/kakorgao/.conda/envs/kat_py_37/

### Read in dynamic inputs from flags using getopts
# Necessary inputs: directory where the hg38 tagaligns are, celltypes list, ABC output dir, exact path to conda env with macs2, number of cores to parallelize
while getopts t:c:o:e:n: flag
do
    case "${flag}" in
        t) tagalign_inprefix=${OPTARG};;
        c) celltypes_fp=${OPTARG};;
        o) ABC_dir=${OPTARG};;
        e) env_fp=${OPTARG};;
        n) num_cores=${OPTARG};;
    esac
done

### Set up steps
# Set up necessary array
mapfile -t celltypes < $celltypes_fp

# Create the hg19 tagAlign directory
tagalign_outdir="${ABC_dir}/inputs/tagAligns_hg19"
if [ ! -d "$tagalign_outdir" ]; then mkdir $tagalign_outdir; fi

# Create the hg19 peaks dir
peaks_outdir="${ABC_dir}/inputs/peak_calls_hg19"
if [ ! -d "$peaks_outdir" ]; then mkdir $peaks_outdir; fi

# Create a log file
log_file="$peaks_outdir/log_file_prepare_peaks_parallel.txt"
if [ ! -f "$log_file" ]; then touch $log_file; fi


### Liftover tagAlign files to hg19
echo "Lifting over tagAligns to hg19 "`date` >> $log_file

# Perform liftover with CrossMap, parallelized with parallel to do all cell types at once
chain="/nfs/lab/ABC/references/hg38ToHg19.over.chain"
parallel -j $num_cores "CrossMap.py bed $chain ${tagalign_inprefix}{}.tagAlign ${tagalign_outdir}/{}.tagAlign" ::: ${celltypes[@]}
echo "" >> $log_file


#### Call peaks on the hg19 tagAlign files
echo "Calling peaks in hg19 "`date` >> $log_file

# Activate necessary conda env
eval "$(conda shell.bash hook)"
conda activate $env_fp
which macs2

# Call peaks with MACS2, parallelized with parallel to do all cell types at once
parallel -j $num_cores "macs2 callpeak -t ${tagalign_outdir}/{}.tagAlign --outdir $peaks_outdir -n {} -q 0.05 --nomodel --keep-dup all -g hs" ::: ${celltypes[@]}


### Gzip all the hg19 tagAlign files to save space
echo "Gzipping tagAlign files "`date` >> $log_file
cd $tagalign_outdir
gzip *
echo "Done "`date` >> $log_file
