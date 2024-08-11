#!/usr/bin/env bash
### In this script we will run the ABC method! This uses the hg19 merged celltype bams, a narrowPeak file (which we will generate here), 
### bulk H3K27ac ChIP-seq data, and ABC's default HiC data (hg19)
### This assumes you have a conda env with the necessary ABC inputs installed. You can copy this one if necessary:
### conda create --name run_abc --clone /home/rlmelton/.conda/envs/final-abc-env

### Read in dynamic inputs from flags using getopts
# Necessary inputs: celltypes list, overall ABC output dir, prefix for output run, filepath to abc env, list of genes "ubiquitously expressed"
# This updated version also needs: prefix for pesudobulk TPM files (for gene exp), and prefix for H3K27ac ChIP-seq bams
# For all prefixes, don't include the / at the end, as this is hardcoded in already
while getopts c:o:r:e:g:x:h: flag
do
    case "${flag}" in
        c) celltypes_fp=${OPTARG};;
        o) ABC_dir=${OPTARG};;
        r) output_prefix=${OPTARG};;
        e) env_fp=${OPTARG};;
        g) gene_list=${OPTARG};;
        x) exp_gene_prefix=${OPTARG};;
        h) chip_prefix=${OPTARG};;
    esac
done


### Set up steps
# Set up necessary array
mapfile -t celltypes < $celltypes_fp

# Create the overall ABC output directory
outdir1="${ABC_dir}/outputs"
if [ ! -d "$outdir1" ]; then mkdir $outdir1; fi
outdir2="${ABC_dir}/outputs/${output_prefix}"
if [ ! -d "$outdir2" ]; then mkdir $outdir2; fi

# Create a log file
log_file="$outdir2/log_file_run_analysis.txt"
if [ ! -f "$log_file" ]; then touch $log_file; fi

# Activate the conda env for running ABC
eval "$(conda shell.bash hook)"
conda activate $env_fp


### Loop through all cell types and run ABC analyses all in a row
for ct in "${celltypes[@]}"
do
  echo "Running ABC scripts on ${ct} cells "`date` >> $log_file
  ct_outdir="${outdir2}/${ct}"
  if [ ! -d "$ct_outdir" ]; then mkdir $ct_outdir; fi

  # Make candidate regions
  echo "Running makeCandidateRegions.py "`date` >> $log_file
  python /path/to/repo/ABC-Enhancer-Gene-Prediction/src/makeCandidateRegions.py \
  --narrowPeak "${ABC_dir}/inputs/peak_calls_hg19/${ct}_peaks.narrowPeak" \
  --bam "${ABC_dir}/inputs/ATAC_bams/hg19_liftover/${ct}.bam" \
  --outDir $ct_outdir \
  --chrom_sizes /path/to/repo/ABC-Enhancer-Gene-Prediction/reference/chr_sizes \
  --regions_blocklist /path/to/repo/ABC-Enhancer-Gene-Prediction/reference/wgEncodeHg19ConsensusSignalArtifactRegions.bed \
  --regions_includelist non-diabetic-islet-multiomics/references/gene_coords.gencodev32.hg19.TSS500bp.bed \
  --peakExtendFromSummit 250 --nStrongestPeaks 150000

  # Run neighborhoods -- add in H3K27ac ChIP-seq and expressed genes (TPM) data here!
  # ALSO CHANGED THE GENES FILE FORM THE PREVIOUS VERISON OF THE SCRIPT
  echo "Running run.neighborhoods.py "`date` >> $log_file
  python /path/to/repo/ABC-Enhancer-Gene-Prediction/src/run.neighborhoods.py \
  --candidate_enhancer_regions "${ct_outdir}/${ct}_peaks.narrowPeak.${ct}.bam.Counts.bed" \
  --genes gene_coords.gencodev32.hg19.TSS500bp.bed/references/gene_coords.gencodev32.hg19.bed \
  --H3K27ac "${chip_prefix}/H3K27ac_${ct}.bam" \
  --ATAC "${ABC_dir}/inputs/ATAC_bams/hg19_liftover/${ct}.bam" \
  --expression_table "${exp_gene_prefix}/${ct}_expressed_genes_values.ALL_GENES.txt" \
  --chrom_sizes /path/to/repo/ABC-Enhancer-Gene-Prediction/reference/chr_sizes \
  --ubiquitously_expressed_genes $gene_list \
  --cellType $ct \
  --outdir "${ct_outdir}/Neighborhoods"

  # Make ABC predictions
  echo "Running predict.py "`date` >> $log_file
  python /path/to/repo/ABC-Enhancer-Gene-Prediction/src/predict.py \
  --enhancers "${ct_outdir}/Neighborhoods/EnhancerList.txt" \
  --genes "${ct_outdir}/Neighborhoods/GeneList.txt" \
  --HiCdir /path/to/downloaded/ABC/BroadHiCdata/ \
  --chrom_sizes /path/to/repo/ABC-Enhancer-Gene-Prediction/reference/chr_sizes \
  --hic_resolution 5000 \
  --scale_hic_using_powerlaw \
  --threshold .02 \
  --cellType $ct \
  --outdir "${ct_outdir}/Prediction" \
  --make_all_putative
  echo "" >> $log_file
done


echo "Done "`date` >> $log_file
