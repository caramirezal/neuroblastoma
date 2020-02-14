## Functions for data processing

liger_save <- function(ligerex,
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


