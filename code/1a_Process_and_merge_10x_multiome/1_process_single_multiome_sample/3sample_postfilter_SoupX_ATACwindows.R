# In this script we will read in the filtered .rds file (from sample_metrics_filtering.R)
# and then perform SoupX on the RNA counts. Afterwards, we will add in the ATAC WINDOWS counts
# (from sample_make_ATAC_lfm_fromBCs.py), regenerate bulk metrics, and perform standard clustering analyses.


### Set up steps
#Read in necessary inputs
args = commandArgs(trailingOnly=TRUE)
sample_dir <- args[1]
output_dir <- args[2] #sample specific subdir of overall outputs dir
reticulate_path <- args[3]
marker_file <- args[4]
short_marker_file <- args[5]
soupx_marker_file <- args[6]

#Set reticulate python
#These are essential reticulate functions that allow us to use python packages in R
#you'll need a conda env with 'leidenalg' and 'pandas' installed to do this
#then route reticulate to the python installed in that conda env with the below functions
Sys.setenv(RETICULATE_PYTHON=sprintf("%s/bin/python",reticulate_path))
library(reticulate)
reticulate::use_python(sprintf("%s/bin/python",reticulate_path))
reticulate::use_condaenv(reticulate_path)

# Load necessary libraries, 
suppressMessages(library(hdf5r))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))
suppressMessages(library(ggpubr))
suppressMessages(library(knitr))
suppressMessages(library(SoupX))

opts_chunk$set(tidy=TRUE)
warnLevel <- getOption('warn')
options(warn = -1)


### Read in .rds Seurat object
print(paste("Reading in the RDS file:",Sys.time()))
adata <- readRDS(file = file.path(output_dir,'intermediate_filtered.rds'))


### Perform SoupX count correction
#Read in raw (tod) and filtered (toc) gex metrics
print(paste("Gathering data and creating the Soup Channel object:",Sys.time()))
wd <- sample_dir
tod <- Seurat::Read10X_h5(file.path(wd, 'raw_feature_bc_matrix.h5'))$`Gene Expression`
DefaultAssay(adata) <- 'RNA'
toc = GetAssayData(object = adata, slot = "counts")

#Compile necessary metadata (pull out the required metadata from the clustered filtered 
#adata object; we need the UMAP coordinates (RD1 and RD2) and the cluster assignments at minimum)
metadata <- (cbind(as.data.frame(adata[["umap.wnn"]]@cell.embeddings),
                   as.data.frame(Idents(adata))))
colnames(metadata) <- c("RD1","RD2","Cluster")

#Make sure the metadata BCs match the toc matrix barcodes
print(paste("Checking all necssary barcode and gene indices match:",Sys.time()))
if (FALSE %in% (row.names(metadata) %in% colnames(toc)) == TRUE) {
  stop('Not all metadata barcodes found in toc')}
if (FALSE %in% (colnames(toc) %in% row.names(metadata)) == TRUE) {
  stop('Not all toc barcodes found in metadata')}

#Make sure the toc and tod indices match
if (FALSE %in% (row.names(tod) %in% row.names(toc)) == TRUE) {
  stop("Not all tod genes found in toc")}
if (FALSE %in% (row.names(toc) %in% row.names(tod)) == TRUE) {
  stop("Not all toc genes found in tod")}
if (FALSE %in% (row.names(toc) == row.names(tod)) == TRUE) {
  stop("toc and tod gene indices don't match")
}

#Create a Soup Channel object and add in the metadata
sc <- SoupChannel(tod,toc)
sc <- setDR(sc,metadata[colnames(sc$toc),c("RD1","RD2")])
sc <- setClusters(sc,setNames(metadata$Cluster,rownames(metadata)))

#Output the estimated percentage of background as a separate file!
background <- sc$soupProfile[order(sc$soupProfile$est,decreasing=TRUE),] 
write.table(background,file.path(output_dir,"SoupX_estimated_background.csv"), sep='\t') 
### keep an eye on this output (some is char, other is #s; not sure if we need stringsAsFactors = T/F?)

#Estimate the contamination fraction with marker gene lists (read these in)
soupx.markers <- read.table(soupx_marker_file, sep='\t', header=TRUE)
Genes1 <- c(soupx.markers$Gene[[1]])
Genes2 <- c(soupx.markers$Gene[[2]])
Genes3 <- c(soupx.markers$Gene[[3]])
Genes4 <- c(soupx.markers$Gene[[4]])

#Select which cells to use with which gene sets to estimate contamination rate
#by adding cluster information we can get better results (applies decisions to
#all cells in a cluster which refines the sets we use)
useToEst = estimateNonExpressingCells(
  sc,
  nonExpressedGeneList = list(A=Genes1, B=Genes2, C=Genes3, D=Genes4))

#Output marker maps of 4 marker gene lists (add titles and use mpar?)
#visualize which cells will be used with which gene sets using marker maps
#blue = expression likely from soup, red = expression likely endogenous
#green outline = will be used with that set, black outline = will not be used
png_path1 <- file.path(output_dir,'plots/SoupX_MarkerMaps.png')
png(png_path1, width = 1152, height = 1152)
theme_set(theme_bw() + 
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                  plot.title = element_text(size = 20, hjust = 0.5),
                  legend.key.size = unit(1, 'cm'), legend.title = element_text(size=15),
                  legend.text = element_text(size=15),
                  axis.title = element_blank(), axis.text = element_text(size=15)))
titles <- paste(soupx.markers$Group, "cell marker list:",soupx.markers$Gene,sep=" ")
plot1 <- plotMarkerMap(sc,geneSet=Genes1,useToEst=useToEst[,"A"]) +
  ggtitle(titles[[1]]) + theme(axis.title.y = element_text(size=15)) 
plot2 <- plotMarkerMap(sc,geneSet=Genes2,useToEst=useToEst[,"B"]) + 
  ggtitle(titles[[2]]) 
plot3 <- plotMarkerMap(sc,geneSet=Genes3,useToEst=useToEst[,"C"]) + 
  ggtitle(titles[[3]]) + theme(axis.title.y = element_text(size=15), axis.title.x = element_text(size=15)) 
plot4 <- plotMarkerMap(sc,geneSet=Genes4,useToEst=useToEst[,"D"]) + 
  ggtitle(titles[[4]]) + theme(axis.title.x = element_text(size=15)) 
figure <- ggarrange(plot1, plot2, plot3, plot4, ncol = 2, nrow = 2,
                    common.legend = TRUE,legend="bottom")
figure
garbage <- dev.off()

#finally use the marker genes and cells to calculate the contamination fraction 
#in each cell and correct the expression profile
sc = calculateContaminationFraction(sc,
                                    list(A=Genes1, B=Genes2, C=Genes3, D=Genes4),
                                    useToEst=useToEst)
print(paste('Estimated global contamination fraction:',100*exp(coef(sc$fit))[[1]],"%"))
out <- adjustCounts(sc)


### Visualize changes to expression profile post SoupX
#Read in short marker gene list and select genes we have data for
marker.genes.short <- scan(short_marker_file,what="",sep="\n")
gene_names1 <- row.names(out)
master_list1 <- marker.genes.short[marker.genes.short %in% gene_names1]
master_list1 <- master_list1[seq(1,9)]

#Output change maps for select marker genes (add titles and use ggarrange)
png_path2 <- file.path(output_dir,'plots/SoupX_ChangeMaps.png')
png(png_path2, width = 1536, height = 1536)
theme_set(theme_bw() + 
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                  plot.title = element_text(size = 20, hjust = 0.5),
                  legend.key.size = unit(1, 'cm'), legend.title = element_text(size=15),
                  legend.text = element_text(size=15),
                  axis.title = element_blank(), axis.text = element_text(size=15)))
plot1 <- plotChangeMap(sc,out,master_list1[[1]]) + ggtitle(master_list1[[1]]) + theme(axis.title.y = element_text(size=15)) 
plot2 <- plotChangeMap(sc,out,master_list1[[2]]) + ggtitle(master_list1[[2]])
plot3 <- plotChangeMap(sc,out,master_list1[[3]]) + ggtitle(master_list1[[3]])
plot4 <- plotChangeMap(sc,out,master_list1[[4]]) + ggtitle(master_list1[[4]]) + theme(axis.title.y = element_text(size=15)) 
plot5 <- plotChangeMap(sc,out,master_list1[[5]]) + ggtitle(master_list1[[5]])
plot6 <- plotChangeMap(sc,out,master_list1[[6]]) + ggtitle(master_list1[[6]])
plot7 <- plotChangeMap(sc,out,master_list1[[7]]) + ggtitle(master_list1[[7]]) + theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))  
plot8 <- plotChangeMap(sc,out,master_list1[[8]]) + ggtitle(master_list1[[8]]) + theme(axis.title.x = element_text(size=15)) 
plot9 <- plotChangeMap(sc,out,master_list1[[9]]) + ggtitle(master_list1[[9]]) + theme(axis.title.x = element_text(size=15)) 
figure <- ggarrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, 
                    ncol = 3, nrow = 3, common.legend = TRUE,legend="bottom")
figure
garbage <- dev.off()


### Load Post SoupX RNA into Seurat Object and add ATAC Windows
print(paste('Creating new Seurat Object with SoupX Counts:', Sys.time()))
adata2 = CreateSeuratObject(out)
adata2[['percent.mt']] <- PercentageFeatureSet(adata2, pattern = '^MT-')

#add in previous raw RNA data as another assay (RNA_raw)
DefaultAssay(adata) <- 'RNA'
raw_rna <-  GetAssayData(object = adata, slot = "counts")
raw_rna_assay <- CreateAssayObject(counts = raw_rna)
adata2[['RNA_raw']] <- raw_rna_assay

#read in the ATAC long format matrix as a table
print(paste("Converting the ATAC long format matrix into a sparse matrix:",Sys.time()))
atac_lfm <- read.table(file.path(output_dir, 'atac.long_fmt.filtered_barcode.mtx.gz'), sep='\t', header=FALSE, stringsAsFactors=FALSE)

#we need to make sure the ATAC BCs match the one currently in adata2$RNA
#double check that all desired BCs are found in the ATAC BCs (and vice versa)
bc_order <- colnames(adata2[["RNA"]])
if (FALSE %in% (bc_order %in% atac_lfm$V2) == TRUE) {
  stop("Not all BCs from adata2 found in ATAC long format matrix")}
if (FALSE %in% (atac_lfm$V2 %in% bc_order) == TRUE) {
  stop("Not all BCs from ATAC long format matrix found in adata2")}

#now set the levels of the lfm based on the desired bc order and reorder based on them
atac_lfm$V2 <- factor(atac_lfm$V2, levels=bc_order)
reordered_lfm <- atac_lfm[order(atac_lfm$V2),]
#do the factoring for windows (V1) in whatever order (automatically alphanumeric?)
reordered_lfm$V1 <- as.factor(reordered_lfm$V1)
#finally make the ATAC sparse matrix and double check the BC order is correct
atac_counts2 <- with(reordered_lfm,sparseMatrix(i=as.numeric(V1), j=as.numeric(V2), x=V3, dimnames=list(levels(V1), levels(V2))))
if (FALSE %in% (colnames(atac_counts2) == bc_order) == TRUE) {
  stop("Unable to create an ATAC sparse matrix with the correct BC order")}

#assemble the other necessary components for a chromatin assay object
print(paste("Creating the Chromatin Assay Object:",Sys.time()))
grange.counts2 <- StringToGRanges(rownames(atac_counts2), sep = c(':', '-'))
grange.use2 <- seqnames(grange.counts2) %in% standardChromosomes(grange.counts2)
atac_counts2 <- atac_counts2[as.vector(grange.use2), ]
suppressMessages(annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86))
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- 'hg38'

#create the chromatin assay object and add to Seurat object
frag.file2 <- file.path(wd, 'atac_fragments.tsv.gz')
chrom_assay2 <- CreateChromatinAssay(counts=atac_counts2, sep=c(':', '-'), genome='hg38', fragments=frag.file2, min.cells=0, min.features=-1, annotation=annotations)
adata2[['ATAC']] <- chrom_assay2


### Compute post-SoupX bulk metrics
#add QC data in
print(paste("Reading in QC data:",Sys.time()))
qc <- read.table(file.path(wd, 'per_barcode_metrics.csv'), sep=',', header=TRUE, stringsAsFactors=1)
qc <- as.data.frame(qc)
rownames(qc) <- qc$gex_barcode
qc <- qc[Cells(adata2), 6:length(colnames(qc))]
adata2 <- AddMetaData(adata2, qc)

#now compute bulk metrics
print(paste("Computing final bulk metrics:",Sys.time()))
print(paste("ATAC median high-quality fragments per cell (atac_fragments):",median(adata2[[]][,'atac_fragments'])))
DefaultAssay(adata2) <- 'ATAC'
adata2 <- TSSEnrichment(adata2)
print(paste("ATAC median TSSe:",median(adata2[[]][,'TSS.enrichment'])))
print(paste("RNA median genes per cell:",median(adata2[[]][,'nFeature_RNA'])))


### Analyze and cluster the post-SoupX object
# RNA analysis
DefaultAssay(adata2) <- 'RNA'
adata2 <- SCTransform(adata2, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims=1:50, reduction.name='umap.rna', reduction.key='rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(adata2) <- 'ATAC'
adata2 <- RunTFIDF(adata2)
adata2 <- FindTopFeatures(adata2, min.cutoff='q0')
adata2 <- RunSVD(adata2)
adata2 <- RunUMAP(adata2, reduction='lsi', dims=2:50, reduction.name='umap.atac', reduction.key='atacUMAP_')

adata2 <- FindMultiModalNeighbors(adata2, reduction.list=list('pca', 'lsi'), dims.list=list(1:50, 2:50))
adata2 <- RunUMAP(adata2, nn.name='weighted.nn', reduction.name='umap.wnn', reduction.key='wnnUMAP_')
adata2 <- FindClusters(adata2, graph.name='wsnn', algorithm=3, resolution = .2, verbose=FALSE)

adata2$log_nCount_ATAC = log(adata2$nCount_ATAC)
adata2$log_nCount_SCT = log(adata2$nCount_SCT)
adata2$log_nFeature_ATAC = log(adata2$nFeature_ATAC)
adata2$log_nFeature_SCT = log(adata2$nFeature_SCT)


### Save final Seurat object to an .rds file 
print(paste("Saving final RDS:",Sys.time()))
saveRDS(adata2, file = file.path(output_dir,'final_filtered.rds'))


### Plotting results (as pngs, more than previous script)
print(paste("Outputting Plots:",Sys.time()))

#3 UMAPs
png_path3 <- file.path(output_dir,'plots/final_filtered_UMAPs.png')
png(png_path3, width = 1728, height = 576)
p1 <- DimPlot(adata2, reduction='umap.rna', group.by='seurat_clusters', label=TRUE, label.size=8, repel=TRUE) + ggtitle('RNA')
p1 <- p1 + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle('RNA only')
p2 <- DimPlot(adata2, reduction='umap.atac', group.by='seurat_clusters', label=TRUE, label.size=8, repel=TRUE) + ggtitle('ATAC')
p2 <- p2 + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle('ATAC only')
p3 <- DimPlot(adata2, reduction='umap.wnn', group.by='seurat_clusters', label=TRUE, label.size=8, repel=TRUE) + ggtitle('WNN')
p3 <- p3 + xlab('UMAP 1')+ ylab('UMAP 2') + ggtitle('Combined')
p1 + p2 + p3 & NoLegend() & theme(plot.title=element_text(hjust=0.5))
garbage <- dev.off()

#Violin plots of metrics
png_path4 <- file.path(output_dir,'plots/final_filtered_metrics.png')
png(png_path4, width = 2304, height = 1536)
p1 <- VlnPlot(adata2, features='nCount_SCT', group.by='seurat_clusters', pt.size=0, log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6, pt.size=0) + geom_hline(yintercept=median(adata2$nCount_SCT), linetype='dashed', lw=2) + 
  theme(plot.title = element_text(size = 25), axis.text = element_text(size=20))
p2 <- VlnPlot(adata2, features='nFeature_SCT', group.by='seurat_clusters', pt.size=0, log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6, pt.size=0) + geom_hline(yintercept=median(adata2$nFeature_SCT), linetype='dashed', lw=2) + 
  theme(plot.title = element_text(size = 25), axis.text = element_text(size=20))
p3 <- VlnPlot(adata2, features='nCount_ATAC', group.by='seurat_clusters', pt.size=0, log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6, pt.size=0) + geom_hline(yintercept=median(adata2$nCount_ATAC), linetype='dashed', lw=2) + 
  theme(plot.title = element_text(size = 25), axis.text = element_text(size=20))
p4 <- VlnPlot(adata2, features='nFeature_ATAC', group.by='seurat_clusters', pt.size=0, log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6, pt.size=0) + geom_hline(yintercept=median(adata2$nFeature_ATAC), linetype='dashed', lw=2) + 
  theme(plot.title = element_text(size = 25), axis.text = element_text(size=20))
p5 <- VlnPlot(adata2, features='nCount_RNA', group.by='seurat_clusters', pt.size=0, log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6, pt.size=0) + geom_hline(yintercept=median(adata2$nCount_SCT), linetype='dashed', lw=2) + 
  theme(plot.title = element_text(size = 25), axis.text = element_text(size=20))
p6 <- VlnPlot(adata2, features='nFeature_RNA', group.by='seurat_clusters', pt.size=0, log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6, pt.size=0) + geom_hline(yintercept=median(adata2$nFeature_SCT), linetype='dashed', lw=2) + 
  theme(plot.title = element_text(size = 25), axis.text = element_text(size=20))
figure <- ggarrange(p1, p3, p5, p2, p4, p6, ncol = 3, nrow = 2)
figure
garbage <- dev.off()

#UMAP feature plots of metrics
png_path5 <- file.path(output_dir,'plots/final_filtered_metrics_UMAP.png')
png(png_path5, width = 1920, height = 480)
p1 <- FeaturePlot(adata2, reduction='umap.wnn', features=c('log_nCount_ATAC','log_nFeature_ATAC','log_nCount_SCT','log_nFeature_SCT'), 
                  cols=c('lightgrey', 'red'), ncol=4) + xlab('UMAP 1') + ylab('UMAP 2')
p1 & theme(plot.title=element_text(hjust=0.5))
garbage <- dev.off()

#Marker gene plot
png_path6 <- file.path(output_dir,'plots/final_filtered_marker_genes.png')
png(png_path6, width = 960, height = 480)
options(repr.plot.width=12, repr.plot.height=6)
marker.genes.long <- scan(marker_file,what="",sep="\n")
p1 <- DotPlot(adata2, assay='SCT', features=marker.genes.long, cluster.idents=TRUE) 
p1 <- p1 + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab('') + ylab('')
p1
garbage <- dev.off()

#Marker gene feature plot -- this currently has workarounds in case no immune cells are found, but not other rare cell types
#Select for short list marker genes that are found in the data
marker.genes.short <- scan(short_marker_file,what="",sep="\n")
gene_names2 <- row.names(adata2[["SCT"]][])
master_list2 <- marker.genes.short[marker.genes.short %in% gene_names2]
master_list2 <- master_list2[seq(1,9)]
sct.master_list2 <- paste('sct_',master_list2,sep="")
print("Marker genes selected for cell type FeaturePlots:")
print(sct.master_list2)

png_path7 <- file.path(output_dir,'plots/final_filtered_marker_genes_UMAP.png')
png(png_path7, width = 1440, height = 1440)
options(repr.plot.width=15, repr.plot.height=15)
p1 <- FeaturePlot(adata2, reduction='umap.wnn', features=sct.master_list2[c(1,2,3)], ncol=3)
p2 <- FeaturePlot(adata2, reduction='umap.wnn', features=sct.master_list2[c(4,5,6)], ncol=3) 
p3 <- FeaturePlot(adata2, reduction='umap.wnn', features=sct.master_list2[c(7,8,9)], ncol=3) 
for (j in c(1,2,3)) {
  p1[[j]] <- p1[[j]] + xlab('')
  p2[[j]] <- p2[[j]] + xlab('')
}
for (j in c(2,3)) {
  p1[[j]] <- p1[[j]] + ylab('') + theme(legend.position = c(.05, .95))
  p2[[j]] <- p2[[j]] + ylab('') + theme(legend.position = c(.05, .95))
  p3[[j]] <- p3[[j]] + ylab('') + theme(legend.position = c(.05, .95))
  p1[[1]] <- p1[[1]] + theme(legend.position = c(.05, .95))
  p2[[1]] <- p2[[1]] + theme(legend.position = c(.05, .95))
  p3[[1]] <- p3[[1]] + theme(legend.position = c(.05, .95))
}
p1 / p2 / p3
garbage <- dev.off()
print(paste("Done:",Sys.time()))
