# Code used to predict cRE-gene links with 3 methods and compare the results
1. Predict cRE-target gene links for a cell type computationally
    1. Run the Single-cell MultimOdal REgulatory Scorer (SMORES) method using paired snRNA-seq and snATAC-seq data: `1a_Run_SMORES.ipynb`
    2. Run the Activity-By-Contact (ABC) model using pseudobulked ATAC-seq data and H3K27ac ChIP-seq data: `1b_ABC_pipeline`
    3. Run the Cicero method using snATAC-seq data: `1c_Run_Cicero_on_Multiome_Object.ipynb`
2. Notebook for processing cRE-gene links from all methods to a common bedpe format: `2_Process_All_Links_to_Common_Format.ipynb`
3. Notebook for combining cRE-gene links from all methods into one file, and also performing basic comparisons between links from different methods: `3_Links_Methods_Basic_Comparisons.ipynb`
4. Notebook for overlapping different subsets of cRE-gene links with eQTLs and calculating the enrichment of concordant eQTLs: `4_Links_Methods_eQTL_Enrichment_Comparisons.ipynb`
5. Notebook for comparing the overlap of cRE-gene links subsets with HiChIP chromosome contact information: `5_Links_Methods_HiChIP_Overlap_Comparisons.ipynb`
6. Notebook for comparing the enrichment of GWAS credible sets in cREs from cRE-gene links: `6_Links_Methods_FINRICH.ipynb`
