# Code used to perform chromatin accessibility QTL analysis 
1. Prepare inputs for RASQUAL
    1. Merge genotypes and visualize genomic PCA analsis for population structure: `00.0_Prepare_vcf_files.ipynb`
    2. Prepare count matrices and generate count PCA's to exclude outliers. Generate covriate matrix: `01_Cell_type_BAM_and_count_matrices.ipynb` 
    3. Add allele specific counts to the vcf using rasqual_tools: `02_prepare_cell_type_vcf_allelic_counts.ipynb`
2. Run RASQUAL:
    1. Run RASQUAL and perform FDR correction using permution background: `03.0_run_caqtls_run1.ipynb`
    2. Subset RASQUAL results by peaks called in the given cell type and redo FDR: `03.1_FIlter_caQTLs_run1.ipynb`
3. Generate genotype class matrix and check enrichment of 2 class vs 3 class variants in significant associations: `04.0_Genotype_class_analysis.ipynb`
4. High level summaries of caQTL results across cell types: `04.1_summaries_run1F.ipynb`
5. Annotate caQTLs based on variants being in the given peak, another, or no peak: `04.2.0_WE_In_Peak.ipynb`
6. Motif analysis:
    1. MotifbreakR analysis to identify enriched of disrupted motifs by cell type: `04.2.1_MotifbreakR_JASPAR.ipynb`
    2. Identification of motifs disrupted by caQTLs with associated variants outside of peaks: `04.2.2_MotifbreakR_JASPAR_Not_In_Peak.ipynb`
7. Compare single cell and bulk:
    1. Concordance with our own "bulk-like" analysis: `04.3_compare_with_Our_Bulk_Run1_ggplot.ipynb`
    2. Concordance with previously published bulk islet caQTLs: `04.3_compare_with_published_islets_caQTLs_Run1_ggplot.ipynb`
8. Concordance and specificity of caQTL results across cell types using mashR: `04.4.1_Interaction_analysis_Run1.ipynb`
9. Added feature specificity to master table of specificity of caQTL features: `04.4.2_Specificity_Table_and_Peaks.ipynb`
10. Run MatrixEQTL for colocalization analysis: `05_WE_100kb_Matrix_QTL_Looped.ipynb`
11. Code for plotting panels associated with figure 3: `06.0_WE_Graphing_Figure_3.ipynb`
