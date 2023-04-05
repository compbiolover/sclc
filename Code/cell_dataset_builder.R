#Name: cell_dataset_builder.R
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Efficiently build cell data set objects
#and get pseudotime information.

cell_dataset_builder <- function(vim_genes = vim,
                                 cell_data = denoised_sc,
                                 cell_meta = gene_metadata,
                                 norm_flag = c("log"),
                                 graph_root = 1.5,
                                 show_traj = TRUE,
                                 point_size = 2.5,
                                 my_root   = "Y_19",
                                 my_monocle_graph = "pseudotime_graph.svg",
                                 my_monocle_graph_genes = "genes_graph.svg",
                                 my_moncole_plot_type = "svg",
                                 my_monocle_plot_dpi = 300,
                                 my_monocle_width = 7.5,
                                 my_monocle_height = 7.5,
                                 my_monocle_unit = "in",
                                 my_cds_filename = "cds_output.rds",
                                 my_pt_data_filename = "pt_data.csv") {
  #Required packages
  require(ggplot2)
  require(monocle3)
  
  #Return list
  cds_return_list <- list()
  
  #Making the CDS object
  vim_genes <- vim_genes
  cell_data_set <- new_cell_data_set(cell_data, gene_metadata = cell_meta)
  cds <- preprocess_cds(cell_data_set,
   num_dim = 100,
    method = "PCA",
    norm_method = norm_flag)
    
  cds <- reduce_dimension(cds, reduction_method = "UMAP")
  cds <- cluster_cells(cds, random_seed = 123)
  partition_plot <- plot_cells(cds, color_cells_by = "partition")
  vim_plot <- plot_cells(cds, genes = vim_genes, cell_size = point_size)
  cds <- learn_graph(cds)
  cds <- order_cells(cds, root_pr_nodes = my_root)
  pt_graph <- plot_cells(cds         = cds,
             color_cells_by         = "pseudotime",
             label_cell_groups      = FALSE,
             label_leaves           = FALSE,
             label_branch_points    = TRUE,
             label_principal_points = TRUE,
             graph_label_size       = graph_root,
             show_trajectory_graph  = show_traj,
             cell_size              = point_size)
  
  #Integrating the new pseudotime data into the count matrix.
  #Here I extract just the pseudotime calculations from the cds
  #object to a use for the swithde function calculation later.
  #The output of this code is a data frame that has 2 columns (pseudotime
  #value and sample name[cell name]) and 375 rows
  #(cells with their associated pseudotime).
  pseudotime_data <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  pseudotime_data <- as.data.frame(pseudotime_data)
  colnames(pseudotime_data) <- "Pseudotime"
  pseudotime_contents <- rownames(pseudotime_data)
  pseudotime_data$Samples <- pseudotime_contents
  
  #Saving the  various plots to .svg or other specified format
  #with particular parameters
  
  #The pseudotime graph
  ggsave(filename = my_monocle_graph,
         plot     = print(pt_graph, newpage = FALSE),
         device   = my_moncole_plot_type, dpi = my_monocle_plot_dpi,
         width    = my_monocle_width, height = my_monocle_height,
         units    = my_monocle_unit)
  
  #The gene expression graph
  ggsave(filename = my_monocle_graph_genes,
         plot     = print(vim_plot, newpage = FALSE),
         device   = my_moncole_plot_type, dpi = my_monocle_plot_dpi,
         width    = my_monocle_width, height = my_monocle_height,
         units    = my_monocle_unit)
  
  #Saving the relevant R objects and data frames
  saveRDS(cds, my_cds_filename)
  write.csv(pseudotime_data, my_pt_data_filename)
  
  #Returning our objects
  cds_return_list[["CDS"]] <- cds
  cds_return_list[["Pseudotime"]] <- pseudotime_data
  cds_return_list[["Cell Progression Graph"]] <- vim_plot
  cds_return_list[["PT Graph"]] <- pt_graph
  return(cds_return_list)
}
