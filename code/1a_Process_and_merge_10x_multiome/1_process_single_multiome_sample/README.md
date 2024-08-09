# 10x Multiome Processing Pipeline
This directory contains scripts and notebooks used to process 10X multiome data of pancreatic islets. The pipeline consists of one bash wrapper script which calls three scripts for processing both the gene expression and chromatin accessibility data. You can find all four scripts in this directory, as well as an explanation of each one and example commands below. **Please note: this pipeline was designed for and run using Seurat v4 and will not work as intended with Seurat v5!**

![Pipeline Diagram](https://github.com/Gaulton-Lab/non-diabetic-islet-multiomics/blob/main/images/General_Multiome_Processing_Pipeline.png)

## Script Information

- [Overall Bash script](https://github.com/Gaulton-Lab/non-diabetic-islet-multiomics/tree/main/code/1a_Process_and_merge_10x_multiome/1_process_single_multiome_sample/multiome_sample_processing.sh): `multiome_sample_processing.sh`
    - Inputs: -s (sample name, used as prefix for outputs), -d (sample data directory, where CellRanger ARC outputs are), -r (path to conda env with reticulate and leidenalg installed), -o (overall output dir)
    - Outputs: log file (`$outs_dir/log_file.txt`)
- [R script 1](https://github.com/Gaulton-Lab/non-diabetic-islet-multiomics/tree/main/code/1a_Process_and_merge_10x_multiome/1_process_single_multiome_sample/1sample_metrics_filtering.R): `1sample_metrics_filtering.R`
    - Inputs: `sample_dir`, `output_dir`, `rna_min_features`, `atac_min_features`, `reticulate_path`, `marker_file`
    - Outputs: `intermediate_filtered.rds`, `filtered_barcodes.txt`, `UMIs_per_BCs.png`, `intermediate_filtered_UMAPs.png`, `intermediate_filtered_metrics.png`, `intermediate_filtered_marker_genes.png`
- [Python script](https://github.com/Gaulton-Lab/non-diabetic-islet-multiomics/tree/main/code/1a_Process_and_merge_10x_multiome/1_process_single_multiome_sample/2ATAC_processing.py): `2ATAC_processing.py`
    - Inputs: —bam (`atac_possorted_bam.bam`), —fragments (`atac_fragments.tsv.gz`),  —keep (`filtered_barcodes.txt`), --blacklist (`/path/to/blacklist.bed`), -outdir
    - Optional inputs (all have default settings): —threads, —memory, —mapquality, —shift, —extsize, —skip-convert, —skip-rmdup, —skip_tagalign, —skip-matrix
    - Outputs: `atac_possorted_bam.compiled.filt.bam`, `atac_possorted_bam.compiled.filt.bam`, `atac_possorted_bam.filt.md.bam`, `atac_possorted_bam.filt.rmdup.bam`, `atac_keep_reads.tagAlign.gz`, `atac_fragments.filtered_barcode.tsv.gz`, `atac_possorted_reads.filtered_barcode.bed.gz`, `atac.long_fmt.filtered_barcode.mtx`
- [R script 2](https://github.com/Gaulton-Lab/non-diabetic-islet-multiomics/tree/main/code/1a_Process_and_merge_10x_multiome/1_process_single_multiome_sample/3sample_postfilter_SoupX_ATACwindows.R): `3sample_postfilter_SoupX_ATACwindows.R`
    - Inputs: `sample_dir`, `output_dir`, `reticulate_path`, `marker_file`, `short_marker_file`, `soupx_marker_file`
    - Outputs: `SoupX_estimated_background.csv`, `SoupX_MarkerMaps.pngs`, `SoupX_ChangeMaps.png`, `final_filtered.rds`, `final_filtered_UMAPS.png`, `final_filtered_metrics.png`, `final_filtered_metrics_UMAP.png`, `final_filtered_marker_genes.png`, `final_filtered_marker_genes_UMAP.png`

### Examples of Running the Bash Script

```bash
sh multiome_sample_processing.sh \
-s R221 \
-d /nfs/lab/projects/multiomic_islet/data/multiomics/cellranger/deep-shallow/R221/outs$
-r /home/hmummey/.conda/envs/reticulate \
-o /nfs/lab/hmummey/multiomic_islet/intermediates/220106_multiome_pipeline_v2_tests
```
