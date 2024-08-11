# Run ABC using 10X snATAC-seq Data

### **DISCLAIMER: this code was developed and used in 2022, and methods to run ABC have majorly changed since then. I'd recommend referring to their [most recent documentation](https://abc-enhancer-gene-prediction.readthedocs.io/en/latest/) instead of using these scripts. We're just including them here for transparency.**

This directory contains scripts and sample files used for running ABC using ATAC-seq data from 10X’s Multiome snATAC-seq. This also incorporates bulk islet H3K27ac ChIP-seq data from [Arda et al. (2018)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79468). 

**NOTES**: 
1. You need to have called peaks by the cell types you want to run ABC on using our internal pipeline before this (or at least have tagAligns split by cell type and a list of all barcodes in your cell type; both are made in the peak calling pipeline). Double check that your input files are named in the correct conventions as well. 
2. The final script which runs ABC also uses some scRNA-seq outputs, a list of genes ubiquitously expressed in your celltypes, and a file with pseudobulked TPM values for all genes per celltype. The files used in this example were created from the snRNA component of data from 10X's Multiome assay. The gene expression input files contribute very little to the ABC runs and could be removed with little impact, if you don't have these inputs, simply modify the `ABC_run_analysis` script to not use the `--ubiquitously_expressed_genes` and/or `--expression_table` flags.
3. This script just needs an input of an ABC run directory and then it will create the appropriate subdirectories to store the processed input files and output files.
4. These scripts are used to run code from the official ABC github repo: [https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction) which we have downloaded to a central location on our lab’s servers that the scripts access. 


### General Workflow

![Pipeline Schematic](https://github.com/Gaulton-Lab/non-diabetic-islet-multiomics/blob/main/images/ABC_pipeline_graphic.png)

![File struction](https://github.com/Gaulton-Lab/non-diabetic-islet-multiomics/blob/main/images/ABC_file_structure.png)

## Step 1 - Process bam files by cell type

To run ABC, we will need pseudobulk ATAC-seq bam files by cell type. To make these we will first separate out the sample bam files by cell type and then merge all sample bams by cell type. Additionally these bam files will be lifted over to hg19 in preparation for running ABC with an hg19 HiC reference. This script uses parallel to process all celltype bams for a sample at once (up to the max number of cores you provide), thus each sample will be processed at the speed it take to process the largest cell type. 

**NOTE:** this step can (and should) be run simultaneously with Step 2, since they are making different input data.

Prepare bams script: `ABC_prepare_bams_parallel.sh`

- Inputs: -s `samples.txt`, -b `/path/to/deduplicated/sample/bams`, -c `celltypes.txt`, -o `/path/to/ABC/run/directory`, -p `/path/to/peak/call/outputs`, -n (number of cores for parallelization)
    - `samples.txt` is a \n delimited .txt file with the names of all samples to use data from
    - `celltypes.txt` is a \n delimited .txt file with the names of all cell types to perform ABC on
- Outputs: a directory of ATAC_bams with subdirectories by sample containing cell type split bam files, a directory of hg19 ATAC bams with merged cell type bam files and subdirectories by sample containing cell type split bam files
    - Ex: `/nfs/lab/projects/islet_ABC/inputs/ATAC/bams`
- Sample command:

```bash
bash ABC_prepare_bams.sh \
-s islet_multiome_samples.txt \
-b /nfs/lab/projects/multiomic_islet/outputs/multiome/dedup_bams \
-c islet_multiome_celltypes.txt \
-o /nfs/lab/projects/islet_ABC \
-p /nfs/lab/projects/multiomic_islet/outputs/multiome/call_peaks/31mergedPeaks_final_majorCTs \
-n 10
```

This script internally calls another R script to perform the initial sectioning of bams by celltypes: `make_celltype_ATAC_bams.R`

- Inputs: directory that contains deduplicated bams (by sample), directory to put cell type bams into, directory of peak calling outputs (containing files with the barcodes of each cell type), `celltypes.txt`, `samples.txt`, number of cores for parallelization
- Outputs: sample bams split by cell type (all in sample subdirectories)

## Step 2 - Call peaks in hg19

Again, because we are running ABC with an hg19 HiC reference, we need to recall peaks in hg19. To do so, we will take the tagAlign files from our hg38 peak calling and lift them over to hg19, then call peaks with MACS2. This script uses parallel to process all celltypes at once (up to the max number of cores you provide).

**NOTE:** this step can (and should) be run simultaneously with Step 1, since they are making different input data.

Prepare peaks script: `ABC_prepare_peaks_parallel.sh`

- Inputs: -t `/path/to/tagAligns/with/fileprefix`, -c `celltypes.txt`, -o `/path/to/ABC/run/directory`, -e `/path/to/conda/env/with/MACS2`, -n (number of cores for parallelization)
    - `/path/to/tagAligns/with/fileprefix` includes everything in the tagAlign file name until `${celltype}.tagAlign` (not gzipped files!)
    - `celltypes.txt` is a \n delimited .txt file with the names of all cell types to perform ABC on
- Outputs: hg19 tagAlign files (gzipped!), and standard MACS2 peak calling outputs (narrowPeak, peaks.xls, summits.bed) for each cell type
    - Ex: `/nfs/lab/projects/islet_ABC/inputs/peak_calls_hg19`
- Sample command:

```bash
bash ABC_prepare_peaks.sh \
-t "/nfs/lab/projects/multiomic_islet/outputs/multiome/call_peaks/31mergedPeaks_final_majorCTs/split." \
-c tests/islet_multiome_celltypes_SMALL.txt \
-o /nfs/lab/projects/islet_ABC \
-e /home/hmummey/.conda/envs/call_peaks \
-n 10
```

## Step 3 - Run ABC

Now that we have all the inputs prepared, we can run ABC! This script just calls the necessary python scripts from the ABC github repo. This script also uses downloaded H3K27ac ChIP-seq data, which is optional when running ABC.

Run ABC script: `ABC_run_analysis_h3k27ac.sh`

- Inputs: -c `celltypes.txt`, -o `/path/to/ABC/run/directory`, -r `prefix_for_run`, -e `/path/to/ABC/conda/env`, -g `expressed_genes_list.txt`, -x `/path/to/TPM/lists`, -h `/path/to/ChIP/hg19/bams`
    - `celltypes.txt` is a \n delimited .txt file with the names of all cell types to perform ABC on
    - `/path/to/tagAligns/with/fileprefix` includes everything in the tagAlign file name until `${celltype}.tagAlign` (not gzipped files!)
    - The outputs will be saved in a subdirectory named after the provided `prefix_for_run`
    - You can create the ABC conda env from the .yaml file in their github repo (`/usr/local/anaconda3/bin/conda env create -f /nfs/lab/ABC/ABC-Enhancer-Gene-Prediction/abcenv.yml`), or just clone it from Rebecca with: `conda create --name run_abc --clone /home/rlmelton/.conda/envs/final-abc-env`
    - The `/path/to/TPM/lists` directory contains files of all genes and their TPM value by celltype in the format: `{celltype}_expressed_genes_values.ALL_GENES.txt`
- Outputs: all ABC outputs and two files with lifted over EnhancerPredictions.txt (`EnhancerPredictions.hg38.txt` and `EnhancerPredictions.genecoords.hg38wgenes.txt`)
    - Ex: `/nfs/lab/projects/islet_ABC/outputs/220809_all_samples2`
- Sample command:

```bash
bash ABC_run_analysis_h3k27ac.sh \
-c tests/islet_multiome_celltypes.txt \
-o /nfs/lab/projects/islet_ABC \
-r "220809_all_samples2" \
-e /home/hmummey/.conda/envs/run_abc \
-g /nfs/lab/projects/multiomic_islet/outputs/multiome/rna_profiles/final_clustering/major_TPM/all_cts_expressed_genes_TPM1.txt \
-x /nfs/lab/projects/islet_ABC/inputs/gex_TPM \
-h /nfs/lab/projects/islet_ABC/inputs/ChIP-seq_bams/hg19_bams 
```

## Step 4 - Liftover the EnhancerPredictions.txt file to hg38

Liftover script: `ABC_hg38_map.sh`

- Inputs: -d `/path/to/directory/with/EnhancerPredictions.txt`
- Outputs: `EnhancerPredictions.hg38.mapped.bedpe`
- Sample command: 

```bash
outdir='/nfs/lab/projects/islet_ABC/outputs/220809_all_samples2/acinar/Prediction'
bash /nfs/lab/ABC/code/ABC_hg38_map.sh -d $outdir
```
