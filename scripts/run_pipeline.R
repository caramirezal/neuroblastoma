## Run pipeline 
## Must be run in the curry cluster
library(liger)
library(Seurat)
library(SeuratWrappers)

## gene expression + atac integration
run_liger <- function(liger,
                      k=10,
                      liger_path,
                      plot_path,
                      interspecies=NULL){
        cat('liger normalization and scaling\n')
        liger <- NormalizeData(liger)
        liger <- FindVariableFeatures(liger)
        liger <- ScaleData(liger, split.by = 'orig.ident', do.center=FALSE)
        
        cat('Running liger\n')
        liger <- RunOptimizeALS(liger, k = k, lambda = 5, split.by = 'orig.ident') 
        liger <- RunQuantileAlignSNF(liger, split.by = 'orig.ident') 
        
        #cat('Run UMAP\n')
        #liger = runUMAP(liger, dims = 1:ncol(liger[['iNMF']]), reduction=='iNMF')
        
        #cat('Plotting\n')
        #pdf(plot_path)
        #DimPlot(liger, group.by = 'orig.ident')
        #dev.off()

        cat('Saving results\n')
        saveRDS(liger, liger_path)
}



cat('Processing liger object\n')
liger <- readRDS('../data/merged_seurat.rds')

cat('Running liger\n')
run_liger(
        liger = liger,
        k = 10,
        liger_path = '../data/batch_correction.rds',
        plot_path = '../figures/batch_correction.pdf',  
        interspecies = TRUE
)

