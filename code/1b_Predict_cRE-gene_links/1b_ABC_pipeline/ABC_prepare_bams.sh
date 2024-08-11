#!/usr/bin/env bash
### In this script we will run necessary preprocessing for ABC, specifically preparing cell type specific bam files
### This will include: separating sample bams by cell type, lifting over to hg19, merging the hg19 bams (and removing the non-merged files)

### Read in dynamic inputs from flags using getopts
# Necessary inputs: sample list, directory where bams are, celltypes list, output dir, directory with peak calling outputs, number of cores
while getopts s:b:c:o:p:n: flag
do
    case "${flag}" in
        s) samples_fp=${OPTARG};;
        b) bam_indir=${OPTARG};;
        c) celltypes_fp=${OPTARG};;
        o) ABC_dir=${OPTARG};;
        p) peak_dir=${OPTARG};;
        n) num_cores=${OPTARG};;
    esac
done


### Set up steps
# Create the overall bam output dir
bam_outdir="${ABC_dir}/inputs/ATAC_bams"
if [ ! -d "$bam_outdir" ]; then mkdir $bam_outdir; fi

# Create a log file
log_file="$bam_outdir/log_file_prepare_bams_parallel.txt"
if [ ! -f "$log_file" ]; then touch $log_file; fi

# Set up necessary arrays
mapfile -t samples < $samples_fp
mapfile -t celltypes < $celltypes_fp

# Get length of samples (for parallelization) and record num_cores as N
num_of_samples=${#samples[@]}
echo "Total number of samples is ${num_of_samples}" >> $log_file
num_of_cts=${#celltypes[@]}
echo "Total number of cell types is ${num_of_cts}" >> $log_file
N=$num_cores


### Make cell type and sample specific bams
make_celltype_bams_Rscript="non-diabetic-islet-multiomics/code/1b_Predict_cRE-gene_links/1b_ABC_pipeline/make_celltype_ATAC_bams.R"
Rscript $make_celltype_bams_Rscript $bam_indir $bam_outdir $peak_dir $celltypes_fp $samples_fp $num_cores >> $log_file
echo "" >> $log_file


### Liftover the cell type and sample specific bams to hg19
echo "Lifting over bams to hg19 "`date` >> $log_file

# Necessary inputs
chain="non-diabetic-islet-multiomics/references/hg38ToHg19.over.chain"
liftover_outdir="${ABC_dir}/inputs/ATAC_bams/hg19_liftover"
if [ ! -d "$liftover_outdir" ]; then mkdir $liftover_outdir; fi

# Loop through every bam file (by sample and celltype) and liftover!
# Uses the parallelize program to run all celltypes for a single sample at once
for s in "${samples[@]}"
do
    echo "Working on sample ${s} "`date` >> $log_file
    out_dir="${liftover_outdir}/${s}" #make output dir for this sample's hg19 bams
    if [ ! -d "$out_dir" ]; then mkdir $out_dir; fi
    parallel -j $num_cores "CrossMap.py bam -a $chain ${ABC_dir}/inputs/ATAC_bams/${s}/${s}_{}.bam ${out_dir}/${s}_{}" ::: ${celltypes[@]}
done
echo "" >> $log_file


### Merge the hg19 bam files by celltype (make sure to NOT put these in a subdir, or will cause issues with using * if you rerun code)
echo "Merging hg19 bams by celltype "`date` >> $log_file

# Merge all celltype bams, sort and index
for (( i=0; i<${num_of_cts}; i++ )); do
    ((j=j%N)); ((j++==0)) && wait
    ct=${celltypes[$i]}
    echo "Working on celltype ${ct} "`date` >> $log_file
    out_fp1="${liftover_outdir}/${ct}.bam"
    all_ct_files="${liftover_outdir}/*/*_${ct}.sorted.bam"
    samtools merge -f -o $out_fp1 $all_ct_files #overwrites any existing file
    out_fp2="${liftover_outdir}/${ct}.sorted.bam"
    sort -k1,1 -k2,2n $out_fp1 > $out_fp2 
    samtools index $out_fp1 &
done


### Remove the hg19_liftover files
cd $liftover_outdir
for s in "${samples[@]}"
do
  rm -r $s*
done
echo "Done "`date`"\n" >> $log_file
