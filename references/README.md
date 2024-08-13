# Reference tables provided
- Basic chromosome information used in processing pipeline: `hg38.chrom.sizes` and `human_autosomes.txt`
- Marker gene lists: `islet_markers_shortlist.txt`, `islet_markers.txt`, and `islet_SoupX_markers.txt`
- Gene coordinates from gencode v32
    - hg19 coordinates: `gene_coords.gencodev32.hg19.bed`
    - hg19 coordinates (only TSS +/- 250bp): `gene_coords.gencodev32.hg19.TSS500bp.bed`
    - hg38 coordinates: `gene_coords.gencodev32.hg38.bed`
    - hg38 coordinates (only TSS +/- 250bp): `gene_coords.gencodev32.hg38.TSS500bp.bed`
- Chain files for liftover with CrossMap: `hg19ToHg38.over.chain` and `hg38ToHg19.over.chain`
- Motif family map (information adapted from TFClass)
    - Full TFClass information table:`220907_WE_Chromvar_to_Gene_Jaspar2022.csv`
    - Family-level summarized table with gene symbols: `220907_WE_Motif_to_Gene_By_Subfam_Complete.JASPAR2022_TFClass.csv`
- Accessible cREs within 100kb of a GWAS lead variant (multiple traits): `final_peaks_mergedPeak.overlap_all_lead_variants_100kb_windows.wINFO.bed`