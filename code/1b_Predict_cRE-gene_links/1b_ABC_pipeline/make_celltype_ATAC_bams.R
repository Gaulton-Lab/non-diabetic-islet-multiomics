# This script is an adaptation of some code Paola sent me to split sample bam files by celltypes. It requires that you have files
# with the barcodes for each celltype (an output of our peak calling pipeline) and is a bit dependent on sample naming conventions.
# It just requires that the sample name and BC are separated by a "_" in the barcodes.txt files.
# This script also assumes you have edit access to the bam files you're going to split, so make sure you do!

# Inputs:
# 1) Directory that contains the bams (preferably deduplicated)
# 2) Directory to put cell type bams into (will create sample subdirs)
# 3) Directory of peak calling outputs (should contain barcodes files)
# 4) File path of text file with all cell types listed (each on new line)
# 5) File path of text file with all samples listed (each on new line)
# 6) Number of cores to use for parallelization

# Outputs:
# 1) Table with barcodes and corresponding celltype, sample name and bam style barcode
# 2) Bam file (and index) for each celltype in each sample (organized into sample subdirs)


### Set up steps
# Read in necessary inputs
args = commandArgs(trailingOnly=TRUE)
bam_indir <- args[1]
bam_outdir <- args[2]
peak_dir <- args[3]
celltypes_fp <- args[4]
samples_fp <- args[5]
num_cores <- args[6]

# Load necessary libraries
suppressMessages(library(stringr))
suppressMessages(library(parallel))


### Create the barcodes dataframe
# Read in barcode lists for each celltype
celltypes   = readLines(celltypes_fp)
barcodelist = data.frame()

# Add necessary info (celltype, sample name and bam style barcode)
for (ct in celltypes){
  bc1          = read.table(file.path(peak_dir,sprintf("%s.barcodes",ct)))
  bc1$celltype = ct
  barcodelist  = rbind(barcodelist, bc1)
}
barcodelist$sample  = stringr::str_match(barcodelist$V1, '(.*)_(.*)')[,2] #assumes sample and barcode are separated by an _, regardless of if sample name has _s
bc                  = str_match(barcodelist$V1, '(.*)_(.*)')
barcodelist$barcode = paste0("CB:Z:", bc[,3])
write.table(barcodelist,file.path(bam_dir,"barcode_celltypes.txt"), quote=F, row.names=F)


### Make separate bam files by sample and celltype (make sure you have edit access to the dedup bam files!)
make_celltype_bam = function(sample, celltypes, bamdir, outdir) {
  print(sprintf('Splitting sample %s bam file: %s', sample, Sys.time()))
  
  # Make an outdir for the sample, save the header, and index the bam
  samp_outdir = file.path(outdir, sample)
  system(paste('mkdir -p', samp_outdir))
  setwd(samp_outdir)
  bam         = file.path(bamdir,sprintf('%s/atac_possorted_bam.filt.rmdup.bam',sample))
  system(paste("samtools view -H", bam ,"> ./SAM_header"))
  system(paste("samtools index", bam ))
  
  for (ct in celltypes) {
    #pull out sample and celltype specific barcodes
    barcodes = barcodelist$barcode[barcodelist$sample==sample & barcodelist$celltype==ct]
    writeLines(barcodes, ct)
    system(paste('samtools view', bam , '| LC_ALL=C grep -F -f', ct , '> filtered_SAM_body'))
    system("cat SAM_header filtered_SAM_body > filtered.sam")
    system(paste0("samtools view -b filtered.sam > ", sample, "_",ct, ".bam"))
    system('rm filtered.sam filtered_SAM_body')
    system(paste0("samtools index ", sample, "_",ct, ".bam"))
  } 
}

# Read in the list of samples and then run make_celltype_bam on each sample in parallel
samples = readLines(samples_fp)
mclapply(samples, function(x) make_celltype_bam(sample=x, celltypes=celltypes, bamdir=bam_indir, outdir=bam_outdir), mc.cores=num_cores)
print(paste("Done with bam celltype splitting: ", Sys.time(), sep=""))