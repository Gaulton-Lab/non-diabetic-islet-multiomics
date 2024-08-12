#!/bin/bash

### Script which prepares an ATAC counts matrix for DESeq for a set of samples and peaks using sinto and featureCounts
### - samples_fp: file with list of samples to run on, each on new line
### - celltypes_fp: file with list of cell types to create featureCount matrices for, each on new line
### - num_cores: number of samples to run in parallel
### - barcode_in_dir: directory with filtered barcode lists by sample (assumes naming convention: sample.filtered_barcodes_wCTs_ofInterest.txt)
### - peaks_saf_fp: file to saf file with peaks to use to make the matrix (made with Mei's script)
### - out_dir: overall directory to write outputs to (script makes subdirs for filtered bams and feature counts outputs)
### - env_with_sinto: path to conda env with sinto installed

### NOTE: because of file naming conventions I had to hardcode in the bam file paths for the HPAP samples... for other uses, modify that one line manually before running this

while getopts s:c:n:b:p:o:e: flag
do
    case "${flag}" in
      s) samples_fp=${OPTARG};;
	  c) celltypes_fp=${OPTARG};;
      n) num_cores=${OPTARG};;
      b) barcode_in_dir=${OPTARG};;
      p) peaks_saf_fp=${OPTARG};;
      o) out_dir=${OPTARG};;
	  e) env_with_sinto=${OPTARG};;
    esac
done

# Read in samples and celltypes
mapfile -t samples < $samples_fp
mapfile -t celltypes < $celltypes_fp

# Set up log file
log_fp=${out_dir}/log.make_featureCounts_mtx.txt

# Filter and divide bam files by cell type 
# First premake the filtered barcode files elsewhere (need second column with cell type assignment)
# Run Sinto in parallel on multiple samples
N=$num_cores
for sample in ${samples[@]}; do
        ((i=i%N)); ((i++==0)) && wait
        (
							echo "Splitting bam file for sample" $sample ":" `date` >> $log_fp
							# HARDCODED BAM FILE PATH -- CHANGE THIS AS NECESSARY
							bam_fp=/nfs/lab/hpap_data/atac/${sample}/Upenn_scATACseq/cellranger_RME/${sample}/outs/possorted_bam.bam
							barcodes_fp=${barcode_in_dir}/${sample}.filtered_barcodes_wCTs_ofInterest.txt
							out_dir=${out_dir}/filt_bams/${sample}
							if [ ! -d "$out_dir" ]; then mkdir -p $out_dir; fi
							${env_with_sinto}/bin/sinto filterbarcodes -b $bam_fp -c $barcodes_fp --outdir $out_dir
							echo "Done splitting bam file for sample" $sample ":" `date` >> $log_fp
        ) &
done
wait

# Prepare inputs for featureCounts
# Make cell type specific bam fp lists
echo "Making cell type-specific bam filepath lists:" `date` >> $log_fp
for celltype in ${celltypes[@]}; do
				cd ${out_dir}/filt_bams
				split_bams_fp=${out_dir}/filt_bams/${celltype}.bam_fps.txt
				touch $split_bams_fp
				for sample in ${samples[@]}; do
								cd $sample
								readlink -f *${celltype}.bam >> $split_bams_fp
								cd ${out_dir}/filt_bams
				done
				#remove any not real file paths from split_bams_fp
				sed -i '/\*/d' $split_bams_fp
done

# Run featureCounts by cell type
fc_dir=${out_dir}/feature_counts
if [ ! -d "$fc_dir" ]; then mkdir -p $fc_dir; fi

for celltype in ${celltypes[@]}; do
				echo "Running featureCounts on" $celltype "cells:" `date` >> $log_fp
				split_bams_fp=${out_dir}/filt_bams/${celltype}.bam_fps.txt
				out_mtx_fp=${fc_dir}/${celltype}_featureCounts_mtx.txt
				bash non-diabetic-islet-multiomics/code/2_Phenotype_associations/2b_quantify_counts_in_peaks/call_featureCounts.sh $split_bams_fp $peaks_saf_fp $out_mtx_fp
				echo "Done running featureCounts on" $celltype "cells:" `date` >> $log_fp
done

echo "All done!" `date` >> $log_fp