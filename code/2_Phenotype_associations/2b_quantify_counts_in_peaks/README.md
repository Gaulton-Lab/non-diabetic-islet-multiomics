# Quantify snATAC-seq counts in a set of peaks

In order to combine the Alberta and HPAP data together for phenotype associations, we need to have the snATAC-seq counts matrices for the exact set of peaks (rows in the matrix). The peaks are from the Alberta object, so here I'll quantify counts for the same set of peaks from the HPAP snATAC-seq data using two packages: `sinto` and `featureCounts`.

## Step 1: Prepare necessary inputs
1. Generate files with the list of barcodes for each sample for cell types of interest where the first column is barcode, second column is cell type with sample name appended to it (Found in notebook `2a_RNA+ATAC_associations_prep.ipynb`)
    1. Format: `${barcode_in_dir}/${sample}.filtered_barcodes_wCTs_ofInterest.txt`
2. Bed file of all peaks you want to quantify counts in

## Step 2: Make a .saf file out of the peaks bed file
Script to convert a bed file to a saf file: `convert_bed_to_saf.py` (written by Mei-lin Okino)
- Inputs: bed file to convert, file path to write saf file to
- Outputs: writes saf file to given file path

Example code:
```
peaks_fp=/nfs/lab/projects/multiomic_islet/outputs/multiome/call_peaks/240304_union_peaks/union_peaks.sort.bed
out_saf=/nfs/lab/projects/multiomic_islet/outputs/multiome/trait_associations/HPAP_ATAC_processing_featureCounts/240308_union_peaks_sort.saf
python convert_bed_to_saf.py $peaks_fp $out_saf
```

## Step 3: Run sinto and featureCounts to quantify counts in peaks
Script to split sample bams by cell type with sinto and then use featureCounts on these to quantify all counts within a set of peaks for a cell type (creates per-donor counts matrix).
- Inputs: 
    - `-s samples_fp`: file with list of samples to run on, each on new line
    - `-c celltypes_fp`: file with list of cell types to create featureCount matrices for, each on new line
    - `-n num_cores`: number of samples to run in parallel
    - `-b barcode_in_dir`: directory with filtered barcode lists by sample (assumes naming convention: sample.filtered_barcodes_wCTs_ofInterest.txt)
    - `-p peaks_saf_fp`: file to saf file with peaks to use to make the matrix (made with Mei's script)
    - `-o out_dir`: overall directory to write outputs to (script makes subdirs for filtered bams and feature counts outputs)
    - `-e env_with_sinto`: path to conda env with sinto installed
- Outputs:
    - `log.make_featureCounts_mtx.txt`: log file
    - `filt_bams`: directory with per donor subdirectories of per-cell type split bam files (from sinto)
    - `feature_counts`: directory with per cell type featureCounts outputs (matrix and summary)

**NOTE: because of file naming conventions I had to hardcode in the bam file paths in this script. Modify this line before running the script: `bam_fp=/nfs/lab/hpap_data/atac/${sample}/Upenn_scATACseq/cellranger_RME/${sample}/outs/possorted_bam.bam`**

Example code:
```
out_dir=/nfs/lab/projects/multiomic_islet/outputs/multiome/trait_associations/240308_HPAP_ATAC_processing/240313_4missing_samples
samples_fp=/nfs/lab/projects/multiomic_islet/outputs/multiome/trait_associations/240308_HPAP_ATAC_processing/hpap_nd_atac_samples_more.txt
celltypes_fp=/nfs/lab/projects/multiomic_islet/outputs/multiome/trait_associations/fin_HPAP_ATAC_processing/celltypes.txt
barcode_in_dir=/nfs/lab/projects/multiomic_islet/outputs/multiome/trait_associations/HPAP_ATAC_processing_featureCounts/sample_barcodes
peaks_saf_fp=/nfs/lab/projects/multiomic_islet/outputs/multiome/trait_associations/240308_HPAP_ATAC_processing/240308_union_peaks_sort.saf
sinto_env=/home/mokino/miniconda3/envs/sinto

bash run_sinto_to_feature_counts.sh -s $samples_fp -c $celltypes_fp -n 6 -b $barcode_in_dir -p $peaks_saf_fp -o $out_dir -e $sinto_env
```
