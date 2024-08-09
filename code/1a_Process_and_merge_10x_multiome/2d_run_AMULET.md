# Instructions for running AMULET on a single 10x Multiome output

1. Run prep script `non-diabetic-islet-multiomics/code/1a_Process_and_merge_10x_multiome/2c_prep_AMULET_inputs.py`
    1. Inputs:
        1. `m`: per_barcode.csv file
        2. `k`: filtered barcodes list output from multiome processing scripts
        3. `o`: file path to save outputs to
    2. Ouput:
        1. `singlecell.csv`
2. Run AMULET (best practice is to run this in a screen)
    1. Sample command:
```bash
/path/to/amulet_zip/AMULET.sh --forcesorted --bambc CB --bcidx 0 --cellidx 8 --iscellidx 9 \
	/path/to/bam.file \
	/path/to/singlecell.csv \
	non-diabetic-islet-multiomics/references/human_autosomes.txt \
	non-diabetic-islet-multiomics/references/blacklist_repeats_segdups_rmsk_hg38.bed \
	/path/to/outputs/ \
	/path/to/amulet_zip/
```
3. Most important output file: `MultipletBarcodes_01.txt` (list of doublet/multiplet barcodes)
