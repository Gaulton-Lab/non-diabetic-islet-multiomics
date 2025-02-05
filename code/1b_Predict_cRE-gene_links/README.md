# Code used to predict cRE-gene links with 3 methods and compare the results
0. Run the Single-cell MultimOdal REgulatory Scorer (SMORES) method using paired snRNA-seq and snATAC-seq data: `0_Run_SMORES.ipynb`
1. Notebook for processing cRE-gene links from all methods to a common bedpe format: `1_Process_All_Links_to_Common_Format.ipynb`
2. Notebook to summarize links from different methods, also creates combined file with all links: `2_Links_Basic_Summary.ipynb`
3. Notebook for performing basic comparisons between links from different methods: `3_Methods_Basic_Comparisons.ipynb`
4. Notebook for overlapping different subsets of cRE-gene links with eQTLs and calculating the enrichment of concordant eQTLs: `4_Method_eQTL_Enrichment_Comparisons.ipynb`
5. Notebook for comparing the overlap of cRE-gene links subsets with HiChIP chromosome contact information: `5_Methods_HiChIP_Overlap_Comparisons.ipynb`
6. Notebook for comparing the enrichment of GWAS credible sets in cREs from cRE-gene links: `6_Methods_FINRICH.ipynb`
