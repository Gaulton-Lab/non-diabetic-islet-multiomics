# Code used to process and merge 10x multiome data
1. Pipeline for processing a single 10x Multiome run, starting with CellRanger ARC count outputs: `1_process_single_multiome_sample`
2. Identify doublets using two methods
    1. Notebook to create inputs for Scrublet from per-sample Seurat objects: `2a_prep_Scrublet_inputs.ipynb`
    2. Notebook to run Scrublet using matrix market files: `2b_run_Scrublet.ipynb`
    3. Python script to convert cellrange arc outputs to necessary input for AMULET: `2c_prep_AMULET_inputs.py`
    4. Example code on how to run AMULET prep script and AMULET itself
3. Notebook for merging all sample pipeline outputs into one Seurat object and removing low quality and doublet barcodes: `3_merge_multiome_samples_and_filter.ipynb`
4. Instructions for calling peaks for each cell type can be found at the following repository: https://github.com/Gaulton-Lab/peak-call-pipeline
5. Notebook for identifying per-donor covariates (biological and technical variables) that are significantly associated with PC representations of the data, and should be accounted for in downstream analyses: `4_find_donor_covars_linked_to_PCs.ipynb`
6. Notebook for using Shannon's entropy to identify genes, cREs and motifs that are specific to one cell type: `5_find_cell_type_specific_comps.ipynb`
