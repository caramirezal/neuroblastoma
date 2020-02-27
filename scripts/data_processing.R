## dependencies
library(liger)
source('scripts/funs.R')
library(Seurat)
library(SeuratWrappers)


## Loading data
file_names <- list.files('data', full.names=TRUE)
f_names <- gsub('data/', '', file_names)
f_names <- gsub('_.*', '', f_names)
f_names <- paste0('nb', f_names)
f_names


## gene expression + atac integration
run_liger <- function(liger,
                      k=10,
                      liger_path,
                      plot_path
                      ){
        cat('Preprocessing data\n')
        liger <- normalize(liger)
        liger <- selectGenes(liger, combine = 'union')
        liger <- scaleNotCenter(liger)
        
        cat('Running liger\n')
        liger <- optimizeALS(liger, k = k) 
        liger <- quantileAlignSNF(liger) 
        
        cat('Run TSNE\n')
        liger = runTSNE(liger)
        
        cat('Plotting visualizations\n')
        pdf(plot_path)
        plotByDatasetAndCluster(liger) #Can also pass in different set of cluster labels to plot
        plotFeature(liger, "nUMI")
        plotWordClouds(liger)
        plotGeneLoadings(liger)
        dev.off()

        cat('Saving results\n')
        saveRDS(liger, liger_path)
}

save_liger <- function(ligerex,
                       dir_path = '.') {
        cat('Creating list object\n')
        liger_res <- list(
                H=ligerex@H,
                cell_data=ligerex@cell.data,
                H_norm=ligerex@H.norm,
                W=ligerex@W,
                V=ligerex@V,
                tsne_coords=ligerex@tsne.coords,
                alignment_clusters=ligerex@alignment.clusters,
                clusters=ligerex@clusters)
        
        cat('Saving results\n')
        saveRDS(liger_res, dir_path)
}

run_liger_seurat <- function(seurat){
     seurat.p <- NormalizeData(seurat)
     seurat.p <- FindVariableFeatures(seurat.p)
     seurat.p <- ScaleData(seurat.p, split.by = "orig.ident", do.center = FALSE)
     seurat.p <- RunOptimizeALS(seurat.p, k = 20, lambda = 5, split.by = "orig.ident")
     seurat.p <- RunQuantileAlignSNF(seurat.p, split.by = "orig.ident")
     seurat.p <- RunUMAP(seurat.p, dims = 1:ncol(seurat.p[["iNMF"]]), reduction = "iNMF")
     saveRDS(seurat.p, 'data/intratumotal_liger_seurat.rds')
}

#############################################################################
##   									   ##
##            Adrenal medula batch correction 				   ##
##									   ##
#############################################################################

#cat('Loading data\n')
#file_names <- list.files('data/adrenal_medulla', full.names=TRUE)
#f_names <- gsub('data/adrenal_medulla/', '', file_names)
#f_names <- gsub('_.*', '', f_names)
#f_names
>>>>>>> d831b7b493c5a2d40aa99bd4891d40077240d509
#counts <- lapply(file_names, readRDS)
#names(counts) <- f_names
#sapply(counts, class)
#sapply(counts, dim)

#cat('Creating liger object\n')
#liger <- createLiger(counts)

#cat('Run liger\n')
#run_liger(
#     liger = liger,
#     k = 10,
#     liger_path = 'data/batch_correction_adrenal_medulla.rds',
#     plot_path = 'figures/batch_correction.pdf'
#)


#liger <- readRDS('data/batch_correction_adrenal_medulla.rds')
#save_liger(liger, 'data/batch_correction_adrenal_medulla_tsne.rds')


########################################################################
##  								      ##
##                 Intratumoral                                       ##
## 								      ##
########################################################################

cat('Reading intratumoral data\n')
seurat <- readRDS('data/intratumoral_merged_seurat.rds')


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

#ligerex <- readRDS('data/batch_correction.rds')

#cat('Plotting visualizations\n')
#pdf('figures/batch_correction.pdf')
#plotByDatasetAndCluster(ligerex) #Can also pass in different set of cluster labels to plot
#plotFeature(ligerex, "nUMI")
#plotWordClouds(ligerex)
#plotGeneLoadings(ligerex)
#dev.off()
=======
cat('Running liger over merged seurta object\n')
run_liger_seurat(seurat)
>>>>>>> d831b7b493c5a2d40aa99bd4891d40077240d509
