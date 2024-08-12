# Code used to perform associations between donor phenotypes and gene expression, cRE accessibility, and motif accessibility
1. Notebook for performing associations between donor covariates and per-donor cell type proportions: `1_Celltype_proportion_associations.ipynb`
2. Prepare necessary inputs for phenotype associations
    1. Notebook for preparing metadata tables, gene counts and cRE accessibility counts tables: `2a_RNA+ATAC_associations_prep.ipynb`
    2. Code use to create per-cell type matrices of snATAC-seq donor counts in a set of peaks (used to get HPAP snATAC-seq counts in Alberta peaks set): `2b_quantify_counts_in_peaks`
3. Notebook for performing RNA associations separately on each dataset using DESeq: `3_RNA_associations_individ.ipynb`
4. Notebook for combining the separate RNA association results with meta-analysis and running pathway analysis on these results: `4_RNA_associations_REM_meta_analysis.ipynb`
5. Perform snATAC-seq phenotype associations
    1. Notebook for performing ATAC associations jointly on both datasets using DESeq: `5a_ATAC_associations_combined.ipynb`
    2. Notebook for performing basic downstream analyses with phenotype-associated cREs: `5b_ATAC_associations_combined_downstream.ipynb`
6. Perform associations between phenotypes and TF motif accessibility
    1. Notebook for running ChromVAR on both Alberta and HPAP datasets separately: `6a_Run_ChromVAR.ipynb`
    2. Notebook for testing for TF motifs associated with phenotypes using a meta analysis of Alberta and HPAP ChromVAR results: `6b_ChromVAR_By_Meta_Interatction_Clean_Looped.ipynb`
