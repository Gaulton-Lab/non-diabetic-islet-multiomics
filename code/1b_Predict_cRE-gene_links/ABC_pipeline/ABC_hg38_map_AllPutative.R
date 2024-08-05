# This is a really basic script that just does some dataframe manipulation to make sure we combine
# the correct rows of the mapped CRE and gene TSS coords files. There is probably a way to do this
# with bash but this was faster for me to implement. This will work with the AllPutative script outputs only!

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
cre_fp = file.path(outdir, 'EnhancerPredictionsAllPutative.crecoords.hg38.txt')
tss_fp = file.path(outdir, 'EnhancerPredictionsAllPutative.genecoords.hg38.txt')
abc_fp = file.path(outdir, 'EnhancerPredictionsAllPutative.rmHeader.txt')

cre_df = read.table(cre_fp, sep='\t')
tss_df = read.table(tss_fp, sep='\t')
abc_df = read.table(abc_fp, sep='\t')

#get row numbers for unmapped rows (check BOTH liftover files)
cre_unmap_fp = file.path(outdir, 'EnhancerPredictionsAllPutative.crecoords.hg38.txt.unmap')
gene_unmap_fp = file.path(outdir, 'EnhancerPredictionsAllPutative.genecoords.hg38.txt.unmap')
unmapped_rows = c()

if (file.exists(cre_unmap_fp)==TRUE & length(readLines(cre_unmap_fp)) > 0){
    cre_unmap_df = read.table(cre_unmap_fp, sep='\t')
    unmapped_rows = c(unmapped_rows,cre_unmap_df$V25)
}
if (file.exists(gene_unmap_fp)==TRUE & length(readLines(gene_unmap_fp)) > 0){
    gene_unmap_df = read.table(gene_unmap_fp, sep='\t')
    unmapped_rows = c(unmapped_rows,gene_unmap_df$V4)
}
unmapped_rows <- unique(unmapped_rows)
#print(unmapped_rows)
print(paste('Number of unmapped rows:',as.character(length(unmapped_rows))))

#make final df by matching using the rowname col previously added
save_rows = abc_df$V25[!abc_df$V25 %in% unmapped_rows]
print(paste('Number of rows in final df:',length(save_rows)))
fin_df = cbind(cre_df[cre_df$V25 %in% save_rows,c(1,2,3)],tss_df[tss_df$V4 %in% save_rows,c(1,2,3)],cre_df[cre_df$V25 %in% save_rows,c(7,21)])
fin_fp = file.path(outdir, 'EnhancerPredictionsAllPutative.hg38.mapped.bedpe')
write.table(fin_df, fin_fp, sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
