library(Seurat)
library(liger)
library(SeuratWrappers)


## Running liger over updated matrix counts of neuroblastoma samples
cat('Loading data\n')
adrenal_seu <- readRDS('data/adrenal_medulla_new/combined_medulla_cc_anno.rds') 
adrenal_ann <- read.table('data/adrenal_medulla_new/cell_type_annotation_medulla_cc.tsv', sep = '\t')
colnames(adrenal_ann) <- gsub('\\.', '', colnames(adrenal_ann))
adrenal_seu <- AddMetaData(adrenal_seu, adrenal_ann)

cat('Running Liger pipeline\n')
ss <- sapply(adrenal_seu$orig.ident, 
             function(i) 
                     ! i %in% c('13667', '13952', '14627', '14742', '14773'))
adrenal_seu.s <- adrenal_seu[,ss]
adrenal_seu.s <- NormalizeData(adrenal_seu.s)
adrenal_seu.s <- FindVariableFeatures(adrenal_seu.s)
adrenal_seu.s <- ScaleData(adrenal_seu.s, split.by = "orig.ident", do.center = FALSE)
adrenal_seu.s <- RunOptimizeALS(adrenal_seu.s, k = 8, lambda = 5, split.by = "orig.ident")
adrenal_seu.s <- RunQuantileAlignSNF(adrenal_seu.s, split.by = "orig.ident")
adrenal_seu.s <- RunUMAP(adrenal_seu.s, dims = 1:ncol(adrenal_seu.s[["iNMF"]]), reduction = "iNMF")
cat('Finish Liger\n')

cat('Saving results\n')
saveRDS(adrenal_seu.s, 'data/liger_ad_med_batch_correction_updated.rds')
