#### get expressed percentage , edited from nichenet ##

library(Seurat)
library(dplyr)
Get.Exp.Genes <- function(seurat_obj, pct=0.1){
 exprs_mat <- seurat_obj@assays$RNA@data
 n_cells_oi_in_matrix <- length(colnames(seurat_obj@assays$RNA))
genes = exprs_mat %>% apply(1, function(x) {
            sum(x > 0)/n_cells_oi_in_matrix
        }) %>% .[. >= pct] %>% names()
 
 
if (n_cells_oi_in_matrix < 5000) {
        genes = exprs_mat %>% apply(1, function(x) {
            sum(x > 0)/n_cells_oi_in_matrix
        }) %>% .[. >= pct] %>% names()
    }
    else {
        splits = split(1:nrow(exprs_mat), ceiling(seq_along(1:nrow(exprs_mat))/100))
        genes = splits %>% lapply(function(genes_indices, exprs, 
            pct, n_cells_oi_in_matrix) {
            begin_i = genes_indices[1]
            end_i = genes_indices[length(genes_indices)]
            exprs = exprs[begin_i:end_i, ]
            genes = exprs %>% apply(1, function(x) {
                sum(x > 0)/n_cells_oi_in_matrix
            }) %>% .[. >= pct] %>% names()
        }, exprs_mat, pct, n_cells_oi_in_matrix) %>% unlist() %>% 
            unname()
}
return(genes)
}
 
