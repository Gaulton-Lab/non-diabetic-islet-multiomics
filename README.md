# Single cell multiome profiling of pancreatic islets reveals physiological changes in cell type-specific regulation associated with diabetes risk.
This repository contains code and reference files used to perform all analyses in our manuscript which describes how both phenotype and genotype change pancreatic islet gene regulation. 

![Project Diagram](https://github.com/Gaulton-Lab/non-diabetic-islet-multiomics/blob/main/images/Project_Diagram.jpeg)

## Data availability
- Check out [multiome.isletgenomics.org](http:multiome.isletgenomics.org) for interactive tools you can use to explore our multimodal single cell map and trait association findings (more tools to come later!)
- Raw data from the paired snRNA-seq and snATAC-seq assays will be available at [GSE273925](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273925) after publication.
- Processed data will be available on the [Common Metabolic Diseases Genome Atlas (CMDGA)](https://cmdga.org/publications/)

## Code
The code contained in this repository is organized as follows:
1. Create a multimodal (snRNA+snATAC) map of gene regulation in islets
    1. Process and clean data from 10x Multiome assays from 28 non-diabetic donors
    2.  Predict cRE-target gene links in six common cell types
3. Test for associations between four phenotypes (age, BMI, HbA1c, and sex) and components of gene regulation
4. Test for associations between genotype and chromatin accessibility (caQTLs)
5. Compare all our findings to known disease risk signals

## Citation
**Please cite any use or adaptation of our analyses at**: [Single cell multiome profiling of pancreatic islets reveals physiological changes in cell type-specific regulation associated with diabetes risk. Mummey, Elison, et al. (2024) bioRxiv](https://www.biorxiv.org/content/10.1101/2024.08.03.606460v1.full#page)
