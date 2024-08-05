#!/usr/bin/env bash

### This is the overall wrapper script for the 10X Multiome data processing pipeline.
### This script was originally designed for the multiomic_islet project and thus has a few references paths
### hardcoded. When applying this script to other data make sure to change these (indicated in line).

### IMPORTANT NOTE: this script relies on the fact that the user has write access to $sample_dir
### make sure you do before running this, otherwise the samtools sort steps in the ATAClfm_Pyscript will fail!


### Read in dynamic inputs from flags using getopts
while getopts s:d:r:o: flag
do
    case "${flag}" in
        s) sample=${OPTARG};;
        d) sample_dir=${OPTARG};;
        r) reticulate_dir=${OPTARG};;
        o) overall_output_dir=${OPTARG};;
    esac
done

### HARDCODED PATHS AND VARIABLES -- CHANGE AS NECESSARY
# locations of the component scripts
filtering_Rscript_fp="/nfs/lab/projects/multiomic_islet/code/indv_sample_processing/1sample_metrics_filtering_v2.R"
ATAClfm_Pyscript_fp="/nfs/lab/projects/multiomic_islet/code/indv_sample_processing/2ATAC_processing_v2.py"
soupX_Rscript_fp="/nfs/lab/projects/multiomic_islet/code/indv_sample_processing/3sample_postfilter_SoupX_ATACwindows_v2.R"

# locations of the marker gene lists
marker_file="/nfs/lab/projects/multiomic_islet/references/islet_markers.txt"
short_marker_file="/nfs/lab/projects/multiomic_islet/references/islet_markers_shortlist.txt"
soupx_marker_file="/nfs/lab/projects/multiomic_islet/references/islet_SoupX_markers.txt"

# BC filtering metrics
min_RNA=500
min_ATAC=1000

# ALSO NOTE THAT THE PYTHON INSTALL TO USE TO RUN THE PYTHON SCRIPT IS CURRENTLY HARDCODED BELOW


### Set up necessary file structure
# make the output directory for this sample
cd $overall_output_dir
outs_dir="$(pwd)/$sample"
if [ ! -d "$outs_dir" ]; then mkdir $outs_dir; fi

# make the plotting subfolder
cd $outs_dir
plots_dir="$(pwd)/plots"
if [ ! -d "$plots_dir" ]; then mkdir $plots_dir; fi

#create a log file
log_file="$outs_dir/log_file.txt"
if [ ! -f "$log_file" ]; then touch $log_file; fi


### Run the scripts in order, with necessary inputs
echo "Individual Sample Processing for Sample "$sample > $log_file
echo "" >> $log_file

# Run the metric filtering script with static thresholds
echo "Starting Rscript #1 "`date` >> $log_file
Rscript $filtering_Rscript_fp $sample_dir $outs_dir $min_RNA $min_ATAC $reticulate_dir $marker_file >> $log_file
echo "" >> $log_file

# Make ATAC long format windows matrix script
cd $sample_dir
echo "Starting Pyscript #1 "`date` >> $log_file
### IMPORTANT NOTE: I've hardcoded in the python install on Ophelia we need to use bc my PATH defaults to
### a different install for some reason. If you're system is setup well you can probably run the commented line instead
#python $ATAClfm_Pyscript_fp --bam "atac_possorted_bam.bam" --keep $outs_dir"/filtered_barcodes.txt" --outdir $outs_dir"/" --fragments "atac_fragments.tsv.gz"  >> $log_file
/usr/bin/python3 $ATAClfm_Pyscript_fp --bam "atac_possorted_bam.bam" --keep $outs_dir"/filtered_barcodes.txt" --outdir $outs_dir --fragments "atac_fragments.tsv.gz"  >> $log_file
echo "" >> $log_file

# Adjust rwx of the 3 output files of ^ (also gzip the lfm matrix file)
chmod -R ugo+rwx "$outs_dir/atac_fragments.filtered_barcode.tsv.gz"
chmod -R ugo+rwx "$outs_dir/atac_possorted_reads.filtered_barcode.bed.gz"
chmod -R ugo+rwx "$outs_dir/atac_possorted_bam.filt.rmdup.bam"
gzip "$outs_dir/atac.long_fmt.filtered_barcode.mtx"
chmod -R ugo+rwx "$outs_dir/atac.long_fmt.filtered_barcode.mtx.gz"

# SoupX and Windows script
cd $outs_dir
echo "Starting Rscript #2 "`date` >> $log_file
Rscript $soupX_Rscript_fp $sample_dir $outs_dir $reticulate_dir $marker_file $short_marker_file $soupx_marker_file >> $log_file
echo "Done "`date` >> $log_file