## dependencies
library(liger)
source('scripts/funs.R')
library(Seurat)



cat('Loading data\n')
file_names <- list.files('data/adrenal_medulla', full.names=TRUE)
f_names <- gsub('data/intratumoral/', '', file_names)
f_names <- gsub('\\.rds', '', f_names)
counts <- lapply(file_names, readRDS)

#cat('Creating Seurat object\n')
#counts_seu <- merge(counts[[1]], counts[2:length(counts)], add.cell.ids = f_names, project = "neuroblastoma")
#merged_dir <- 'data/merged_seurat.rds'
#cat('Saving results to', merged_dir,'\n')
#saveRDS(counts_seu, merged_dir)

merged_seu <- readRDS('data/merged_seurat.rds')
merged_seu
cat('Creating liger object \n')
ligerex <- NormalizeData(merged_seu)

#cat('Preprocessing \n')
#ligerex <- normalize(ligerex)
#ligerex <- selectGenes(ligerex, combine = 'union')
#ligerex <- scaleNotCenter(ligerex)

#cat('Running liger\n')
#ligerex <- optimizeALS(ligerex, k = 10) 
#ligerex <- quantileAlignSNF(ligerex) 
                
#cat('Run TSNE\n')
#ligerex = runTSNE(ligerex)

#cat('Saving results\n')
#saveRDS(ligerex, 'data/batch_correction.rds')

#ligerex <- readRDS('data/batch_correction.rds')

#cat('Plotting visualizations\n')
#pdf('figures/batch_correction.pdf')
#plotByDatasetAndCluster(ligerex) #Can also pass in different set of cluster labels to plot
#plotFeature(ligerex, "nUMI")
#plotWordClouds(ligerex)
#plotGeneLoadings(ligerex)
#dev.off()

#liger_save(ligerex, 'data/batch_correction_sum.rds')
