library(liger)

## Loading data
#file_names <- list.files('data', full.names=TRUE)
#f_names <- gsub('data/', '', file_names)
#f_names <- gsub('_.*', '', f_names)
#f_names <- paste0('nb', f_names)
#f_names

#counts <- lapply(file_names, readRDS)
#lapply(counts, function(c) c[1:1, 1:2] )
#names(counts) <- f_names

#cat('Creating liger object \n')
#ligerex <-  createLiger(counts)

#cat('Preprocessing \n')
#ligerex <- normalize(ligerex)
#ligerex <- selectGenes(ligerex, combine = 'union')
#ligerex <- scaleNotCenter(ligerex)

#cat('Running liger\n')
#ligerex <- optimizeALS(ligerex, k = 10) 
#ligerex <- quantileAlignSNF(ligerex) 
                
cat('Run TSNE\n')
#ligerex = runTSNE(ligerex)

cat('Saving results\n')
#saveRDS(ligerex, 'data/batch_correction.rds')

ligerex <- readRDS('data/batch_correction.rds')

cat('Plotting visualizations\n')
pdf('figures/batch_correction.pdf')
plotByDatasetAndCluster(ligerex) #Can also pass in different set of cluster labels to plot
plotFeature(ligerex, "nUMI")
plotWordClouds(ligerex)
plotGeneLoadings(ligerex)
dev.off()
