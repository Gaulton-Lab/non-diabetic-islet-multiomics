#!/usr/bin/env python3

### This script prepares the necessary input files for AMULET from the CellRanger ARC count outputs
### Code initially developed by Katha Korgaonkar, reorganized by Hannah Mummey

import pandas 
import argparse

def make_amulet_input_csv(args):
    #reading in the per_barcode_metrics file from cellranger arc outputs
    metrics = pandas.read_csv(args.metrics, sep=',')

    #dropping columns we don't need from cellranger arc output
    metrics = metrics.drop(columns= ['gex_barcode', 'atac_barcode', 'excluded_reason',
        'gex_raw_reads', 'gex_mapped_reads', 'gex_conf_intergenic_reads',
        'gex_conf_exonic_reads',  'gex_conf_intronic_reads',
        'gex_conf_exonic_unique_reads', 'gex_conf_exonic_antisense_reads',
        'gex_conf_exonic_dup_reads', 'gex_exonic_umis',
        'gex_conf_intronic_unique_reads', 'gex_conf_intronic_antisense_reads',
        'gex_conf_intronic_dup_reads', 'gex_intronic_umis',
        'gex_conf_txomic_unique_reads', 'gex_umis_count', 'gex_genes_count'])

    #adding in columns found in cellranger atac but not cellranger arc outputs
    metrics['DNase_sensitive_region_fragments'] = 0
    metrics['enhancer_region_fragments'] = 0
    metrics['promoter_region_fragments'] = 0
    metrics['on_target_fragments'] = metrics['atac_TSS_fragments']
    metrics['blacklist_region_fragments'] = 0

    #reading in barcodes from our own "is cell" calls
    keep = pandas.read_csv(args.keep, header = None)
    keep.head()

    #creating and setting new is_cell column
    new_is_cell = []
    for i in metrics['barcode']:
        if i in keep[0].values:
            new_is_cell.append(1)
        elif i not in keep[0].values:
            new_is_cell.append(0)
    metrics['is_cell'] = new_is_cell

    #replacing cellid column
    numid = 1
    cellid = []
    for i in metrics['is_cell']:
        if i == 0:
            cellid = cellid + ["None"]
        elif i == 1:
            cellid = cellid + ["_cell_" + str(numid)]
            numid = numid + 1    
    metrics['cellid'] = cellid

    #reordering columns to matching cellranger atac output (single_cell.csv)
    order = ['barcode','atac_raw_reads','atac_dup_reads', 'atac_chimeric_reads', 'atac_unmapped_reads','atac_lowmapq',
            'atac_mitochondrial_reads','atac_fragments', 'cellid',  'is_cell','atac_TSS_fragments',
            'DNase_sensitive_region_fragments', 'enhancer_region_fragments','promoter_region_fragments', 
            'on_target_fragments','blacklist_region_fragments','atac_peak_region_fragments','atac_peak_region_cutsites']
    metrics = metrics[order]

    #save output
    metrics.to_csv(args.outfp,index=False)


def main(args):
    make_amulet_input_csv(args)
    return


def process_args():
    parser = argparse.ArgumentParser(description='Make AMULET input')
    gen_group = parser.add_argument_group('General arguments')
    gen_group.add_argument('-m', '--metrics', required=True, type=str, help='Path to per_barcode_metrics.csv file')
    gen_group.add_argument('-k', '--keep', required=True, type=str, help='List of barcodes to keep')
    gen_group.add_argument('-o', '--outfp', required=True, type=str, help='File path to save amulet results to')


if __name__ == '__main__':
    args = process_args()
    main(args)