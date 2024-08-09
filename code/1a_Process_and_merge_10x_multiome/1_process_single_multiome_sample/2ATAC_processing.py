#!/usr/bin/env python3

### In this script I perform perform a variety of important processing steps for the multiome ATAC data:
### First I clean up the raw bam files and create deduplicated bams (necessary for QTL analyses).
### Then I create tagAlign files from the raw bam file for a specified list of barcodes (necessary for peak calling).
### Finally, I create a genome wide windows set and quantify the ATAC counts per barcode (also from the specified
### list) in each window and create a long format matrix (used during clustering).


import os
import gzip
import argparse
import subprocess
import pysam
import numpy as np
import pandas as pd
import scipy.sparse
from multiprocessing import Pool
from datetime import datetime


def convert_10X_bam(args):
    ### Convert the 10X bams so that the query names include the barcodes
    print('Converting the bam file:',datetime.now())
    bf = pysam.AlignmentFile(args.bam, 'rb')
    cf = pysam.AlignmentFile(args.outdir + '/atac_possorted_bam.compiled.filt.bam', 'wb', template=bf)

    #Read through raw bam file and add the barcode to the query name
    for read in bf:
        try:
            barcode = read.get_tag('CB').split('-1')[0]
        except:
            pass
            continue
        read.query_name = barcode + '_' + read.query_name
        cf.write(read)
    bf.close()
    cf.close()
    return


def remove_duplicate_reads(args):
    ### Remove duplicate reads from bam files
    print('Removing duplicate reads:',datetime.now())
    filt_bam = args.outdir + '/atac_possorted_bam.compiled.filt.bam'
    markdup_bam = args.outdir + '/atac_possorted_bam.filt.md.bam'
    rmdup_bam = args.outdir + '/atac_possorted_bam.filt.rmdup.bam'

    #Set up methods to process the bam files to remove duplicates
    filt_cmd = ['samtools', 'view', '-bu', '-q', str(args.mapquality), '-F', '256', '-F', '512', '-F', '2048', filt_bam] #filter bams
    sortname_cmd = ['samtools', 'sort', '-n', '-m' , str(args.memory)+'G', '-@', str(args.threads), '-'] #sort bams by read names
    fixmate_cmd = ['samtools', 'fixmate', '-r', '-', '-'] #fill in mate coordinates, remove secondary and unmapped reads
    sortpos_cmd = ['samtools', 'sort', '-m', str(args.memory)+'G', '-@', str(args.threads), '-o', markdup_bam] #sort bams by leftmost coordinates
    index_cmd = ['samtools', 'index', markdup_bam] #index ^ file
    rmdup_cmd = ['samtools', 'view', '-@', str(args.threads), '-b', '-f', '3', '-F', '1024', markdup_bam] #filter bams for those mapped with with pair, no PCR duplicates
    rmdup_cmd.extend(['chr{}'.format(c) for c in list(map(str, range(1,23))) + ['X','Y']]) #add chr prefix

    #Execute processing methods
    with open(os.devnull, 'w') as null:
        filt = subprocess.Popen(filt_cmd, stdout=subprocess.PIPE)
        sortname = subprocess.Popen(sortname_cmd, stdin=filt.stdout, stdout=subprocess.PIPE)
        fixmate = subprocess.Popen(fixmate_cmd, stdin=sortname.stdout, stdout=subprocess.PIPE)
        subprocess.call(sortpos_cmd, stdin=fixmate.stdout)
        subprocess.call(index_cmd)
    if os.path.isfile(markdup_bam): #only execute this if previous steps were successful in creating the markdup_bam file
        with open(rmdup_bam, 'w') as bam_out:
            subprocess.call(rmdup_cmd, stdout=bam_out)
    else:
        raise FileNotFoundError('{} not found!'.format(markdup_bam))
    return


def create_tagAlign(args):
    ### Filter the bam file for keep barcodes and output in tagAlign format (necessary for peak calling)
    print('Creating the tagAlign file:',datetime.now())

    #Read in the keep barcodes and store them as a set
    keep = open(args.keep).read().splitlines()
    keep_set = set(keep)
    tagalign_file = args.outdir + '/atac_keep_reads.tagAlign.gz'

    #Read through the raw bamfile and add lines to tagalign file if passes filters and is in keep_set barcodes
    bamfile = pysam.AlignmentFile(args.bam, 'rb')
    with gzip.open(tagalign_file, 'wt') as tagout:
        genome_size = {item['SN']:item['LN'] for item in bamfile.header['SQ']}
        #Apply filters to each read from the raw bamfile
        for read in bamfile:
            if not read.is_proper_pair or read.is_duplicate or read.is_secondary or read.is_qcfail or read.is_supplementary:
                continue
            if read.mapping_quality < args.mapquality:
                continue
            if not read.has_tag('CB'):
                continue
            #Get barcode and see if it's in the keep set, if not don't proceed with this read
            barcode = read.get_tag('CB')
            if not barcode in keep_set:
                continue
            read_chr = read.reference_name
            if not read_chr in ['chr{}'.format(i) for i in range(1,23)] + ['chrX','chrY','chrM']:
                continue
            #Shift and extend the read position
            read_start = max(1, read.reference_end - args.shift - args.extsize - 5 if read.is_reverse else read.reference_start + args.shift + 4)
            read_end = min(genome_size[read_chr], read.reference_end - args.shift - 5 if read.is_reverse else read.reference_start + args.shift + args.extsize + 4)
            read_qual = read.mapping_quality
            if read.is_reverse:
                read_orient = '-'
            else:
                read_orient = '+'
            #Add the read to tagAlign file in correct format
            print(read_chr, read_start, read_end, barcode, read_qual, read_orient, sep='\t', file=tagout)
        bamfile.close()
        return


def create_ATAC_windows_lfm(args):
    ### Use the keep barcodes list and a genome wide windows list to generate a
    ### long format matrix (lfm) of windows x barcodes for the ATAC data

    # Read in the barcodes keep list and make a filtered fragments file
    print("Creating the filtered fragment file:",datetime.now())
    keep_set = open(args.keep).read().splitlines()
    frags = pd.read_table(args.fragments, sep='\t', header=None , skiprows=51)
    frags = frags.loc[frags[3].isin(keep_set)]
    frags.to_csv(args.outdir+'/atac_fragments.filtered_barcode.tsv.gz', sep='\t', header=False, index=False, compression='gzip')

    # Go through the raw bam file and output as a filtered bed file (filtered on quality only for now)
    print("Creating the filtered bed file:",datetime.now())
    reads_file = args.outdir + '/atac_possorted_reads.filtered_barcode.bed.gz'
    if not os.path.isfile(reads_file):
        bamfile = pysam.AlignmentFile(args.bam, 'rb')
        with gzip.open(reads_file, 'wt') as bedout:
            genome_size = {item['SN']:item['LN'] for item in bamfile.header['SQ']}
            for read in bamfile:
                if not read.is_proper_pair or read.is_duplicate or read.is_secondary or read.is_qcfail or read.is_supplementary:
                    continue
                if read.mapping_quality < 30:
                    continue
                if not read.has_tag('CB'):
                    continue
                read_chr = read.reference_name
                read_start = max(1, read.reference_end - args.shift - args.extsize - 5 if read.is_reverse else read.reference_start + args.shift + 4)
                read_end = min(genome_size[read_chr], read.reference_end - args.shift - 5 if read.is_reverse else read.reference_start + args.shift + args.extsize + 4)
                read_barcode = read.get_tag('CB')
                if not read_chr in ['chr{}'.format(i) for i in range(1,23)] + ['chrX','chrY','chrM']:
                    continue
                print(read_chr, read_start, read_end, read_barcode, sep='\t', file=bedout)

    # Sort the bed file and only keep reads with barcodes in the keep set (modifies file in place)
    print("Sorting the filtered bed file:",datetime.now())
    reads = pd.read_table(reads_file, sep='\t', header=None)
    reads = reads.sort_values([0,1]).loc[reads[3].isin(keep_set)]
    reads = reads.loc[reads[0]!='chrM']
    reads.to_csv(reads_file, sep='\t', header=False, index=False)

    # Use bed tools to generate a windows list and intersect it with the sorted and filtered bed file, then print out as a long format matrix
    print("Creating the windows LFM:",datetime.now())
    mkw = subprocess.Popen(['bedtools', 'makewindows', '-g', '/nfs/lab/elisha/nPOD_output/scripts/references/hg38.chrom.sizes', '-w', '5000'], stdout=subprocess.PIPE) #make 5kb windows from the genome
    bl = subprocess.Popen(['bedtools', 'intersect', '-a' , '-', '-b', '/nfs/lab/ref/hg38-blacklist.v3.bed', '-v'], stdin=mkw.stdout, stdout=subprocess.PIPE) #remove windows that are in the blacklist
    pair = subprocess.Popen(['bedtools', 'intersect', '-a', reads_file, '-b', '-', '-wa', '-wb'], stdin=bl.stdout, stdout=subprocess.PIPE) #intersect the reads file with the windows, keep those that overlap a window
    awk = subprocess.Popen(['awk', 'BEGIN{{FS=OFS=\"\\t\"}} {{print $5\"-\"$6\"-\"$7,$4}}'], stdin=pair.stdout, stdout=subprocess.PIPE) #reformat the window name
    sort = subprocess.Popen(['sort', '-S', '48G'], stdin=awk.stdout, stdout=subprocess.PIPE) #sort the reads by window they overlap
    uniq = subprocess.Popen(['uniq', '-c'], stdin=sort.stdout, stdout=subprocess.PIPE) #count the number of reads that overlap each window for each BC
    with open(args.outdir + '/atac.long_fmt.filtered_barcode.mtx', 'w') as f_out:
        #print('peak', 'barcode', 'value', sep='\t', file=f_out)
        subprocess.call(['awk', 'BEGIN{{OFS=\"\\t\"}} {{print $2,$3,$1}}'], stdin=uniq.stdout, stdout=f_out)
    return


def main(args):
    if not args.skip_convert:
        convert_10X_bam(args)
    if not args.skip_rmdup:
        remove_duplicate_reads(args)
    if not args.skip_tagalign:
        create_tagAlign(args)
    if not args.skip_matrix:
        create_ATAC_windows_lfm(args)
    print("Done:",datetime.now())
    return


def process_args():
    parser = argparse.ArgumentParser(description='Filter fragments file')
    io_group = parser.add_argument_group('I/O arguments')
    io_group.add_argument('-b', '--bam', required=True, type=str, default='atac_possorted_bam.bam', help='Path to 10X output ATAC bam file')
    io_group.add_argument('-k', '--keep', required=True, type=str, help='List of barcodes to keep')
    io_group.add_argument('-o', '--outdir', required=True, type=str, help='Desired output directory')
    io_group.add_argument('-f', '--fragments', required=False, type=str, default='atac_fragments.tsv.gz', help='Path to 10X output ATAC fragments file')

    align_group = parser.add_argument_group('Alignment arguments')
    align_group.add_argument('-t', '--threads', required=False, type=int, default=8, help='Number of threads to use for alignment [8]')
    align_group.add_argument('-m', '--memory', required=False, type=int, default=4, help='Maximum amount of memory (G) per thread for samtools sort [4]')
    align_group.add_argument('-q', '--mapquality', required=False, type=int, default=30, help='Mapping quality score filter for samtools [30]')

    tagalign_group = parser.add_argument_group('tagAlign generation arguments')
    tagalign_group.add_argument('-s', '--shift', required=False, type=int, default=-100, help='Read shift length')
    tagalign_group.add_argument('-e', '--extsize', required=False, type=int, default=200, help='Read extension size')

    skip_group = parser.add_argument_group('Skip steps')
    skip_group.add_argument('--skip-convert', required=False, action='store_true', default=False, help='Skip bam conversion step')
    skip_group.add_argument('--skip-rmdup', required=False, action='store_true', default=False, help='Skip duplicate removal step')
    skip_group.add_argument('--skip-tagalign', required=False, action='store_true', default=False, help='Skip generate tagAlign step')
    skip_group.add_argument('--skip-matrix', required=False, action='store_true', default=False, help='Skip long format matrix generation step')
    return parser.parse_args()


if __name__ == '__main__':
    args = process_args()
    main(args)
