# This is a really basic script that just does some dataframe manipulation to make sure we combine
# the correct rows of the mapped CRE and gene TSS coords files.

# Inputs:
# 1) ABC Prediction directory for one celltype run (made by our pipeline and the ABC scripts)

# Outputs:
# 1) Bedpe file of hg38 mapped CRE coords, gene TSS coords, gene name, ABC score (unmapped rows removed)


### Set up steps
# Read in necessary inputs
args = commandArgs(trailingOnly=TRUE)
outdir <- args[1]

### Combine the relevant dataframes with the correct rows
#read in all mapped (and original unmapped) files
cre_fp = file.path(outdir, 'EnhancerPredictions.crecoords.hg38.txt')
tss_fp = file.path(outdir, 'EnhancerPredictions.genecoords.hg38.txt')
abc_fp = file.path(outdir, 'EnhancerPredictions.rmHeader.txt')

cre_df = read.table(cre_fp, sep='\t')
tss_df = read.table(tss_fp, sep='\t')
abc_df = read.table(abc_fp, sep='\t')

#get row numbers for unmapped rows (check BOTH liftover files)
cre_unmap_fp = file.path(outdir, 'EnhancerPredictions.crecoords.hg38.txt.unmap')
gene_unmap_fp = file.path(outdir, 'EnhancerPredictions.genecoords.hg38.txt.unmap')
unmapped_rows = c()

if (file.exists(cre_unmap_fp)==TRUE & length(readLines(cre_unmap_fp)) > 0){
    cre_unmap_df = read.table(cre_unmap_fp, sep='\t')
    unmapped_rows = c(unmapped_rows,cre_unmap_df$V9)
}
if (file.exists(gene_unmap_fp)==TRUE & length(readLines(gene_unmap_fp)) > 0){
    gene_unmap_df = read.table(gene_unmap_fp, sep='\t')
    unmapped_rows = c(unmapped_rows,gene_unmap_df$V4)
}
unmapped_rows <- unique(unmapped_rows)
#print(unmapped_rows)
print(paste('Number of unmapped rows:',as.character(length(unmapped_rows))))

#make final df by matching using the rowname col previously added
save_rows = abc_df$V9[!abc_df$V9 %in% unmapped_rows]
print(paste('Number of rows in final df:',length(save_rows)))
fin_df = cbind(cre_df[cre_df$V9 %in% save_rows,c(1,2,3)],tss_df[tss_df$V4 %in% save_rows,c(1,2,3)],cre_df[cre_df$V9 %in% save_rows,c(5,8)])
fin_fp = file.path(outdir, 'EnhancerPredictions.hg38.mapped.bedpe')
write.table(fin_df, fin_fp, sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
