# Name: cell_dataset_builder.R
# Author: Andrew Willems <awillems@vols.utk.edu>
# Purpose: Efficiently build cell data set objects and get pseudotime information.
# Documentation:

# Required packages
# ggplot2: for plotting
# igraph: for building graphs
# monocle3: for building cell data set objects and ordering cells by pseudotime

# Function arguments
# The cell_dataset_builder function takes several arguments:
#   
# gene_expression: a character vector specifying the gene expression to plot.
# gene_title: a character string specifying the gene expression to use in plot titles.
# cell_data: a matrix of gene expression data, with rows as cells and columns as genes.
# gene_meta: a dataframe with gene metadata.
# cell_meta: a dataframe with cell metadata.
# norm_flag: a character string specifying the normalization method to use.
# umap_seed: an integer specifying the seed for the UMAP algorithm.
# cluster_res: a numeric value specifying the resolution for clustering cells.
# n_dims: an integer specifying the number of dimensions to use for dimension reduction (PCA only).
# reduction_method: a character string specifying the dimension reduction method to use.
# reduce_prep_method: a character string specifying the dimension reduction method for the preprocessing method.
# preprocess_method: a character string specifying the preprocessing method to use after dimension reduction.
# graph_root: a numeric value specifying the position of the root node in the pseudotime ordering graph.
# show_traj: a logical value indicating whether to show the trajectory graph.
# point_size: a numeric value specifying the size of the data points in the plots.
# pt_root: a character string specifying the name of the root node in the pseudotime ordering graph.
# pt_title: a character string specifying the pseudotime label to use in plot titles.
# plot_file: a character string specifying the filename to save the expression graph plot to.
# plot_type: a character string specifying the type of file to save the plot as.
# use_ordering_funct: a logical value indicating whether to use a function to order cells by pseudotime.
# plot_dpi: an integer specifying the resolution to save the plot as.
# plot_width: a numeric value specifying the width of the plot.
# plot_height: a numeric value specifying the height of the plot.
# plot_unit: a character string specifying the unit of measurement for the plot dimensions.
# cds_filename: a character string specifying the filename to save the CDS object as.
# pt_file: a character string specifying the filename to save the pseudotime plot to.
# pt_data_filename: a character string specifying the filename to save the pseudotime data to.
 
# Function output
# The cell_dataset_builder function returns a list with the following elements:
#   
# CDS: the cell data set object.
# Cell_Progression_Graph: the expression graph plot.
# PT_Graph: the pseudotime plot.


# Required packages
library(ggplot2)
library(igraph)
library(monocle3)

cell_dataset_builder <- function(
    gene_expression = vim,
    gene_title = "VIM",
    cell_data = denoised_sc,
    gene_meta = gene_metadata,
    cell_meta = cell_metadata,
    norm_flag = "none",
    umap_seed = 123,
    cluster_res = 1e-5,
    n_dims = 100,
    reduction_method = "UMAP",
    reduce_prep_method = "PCA",
    preprocess_method = "PCA",
    graph_root = 1.5,
    show_traj = TRUE,
    point_size = 2.5,
    pt_root = "Y_19",
    pt_title = "TUBA1A",
    plot_file = "expression_graph.png",
    plot_type = "png",
    use_ordering_funct = TRUE,
    plot_dpi = 600,
    plot_width = 7,
    plot_height = 7,
    plot_unit = "in",
    cds_filename = "cds.rds",
    pt_file = "pt_graph.png",
    pt_data_filename = "pt_data.csv") {
  
  # Make the CDS object
  cds <- new_cell_data_set(
    expression_data = cell_data,
    gene_metadata = gene_meta,
    cell_metadata = cell_meta
  )

  cds <- preprocess_cds(cds,
    num_dim = n_dims, method = preprocess_method,
    norm_method = norm_flag
  )
  
  cds <- reduce_dimension(cds,
    reduction_method = reduction_method,
    preprocess_method = reduce_prep_method
  )
  plot_cells(cds, genes = gene_expression)
  cds <- cluster_cells(cds, random_seed = umap_seed, resolution = cluster_res)
  cds <- learn_graph(cds)


  # Function to calculate the lowest expression point for the gene and select
  # it as the root of pseudotime ordering
  get_earliest_principal_node <- function(cds) {
    cell_ids <- which.min(colData(cds_test)[, "Size_Factor"])

    closest_vertex <-
      cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <-
      igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
      (which.max(table(closest_vertex[cell_ids, ]))))]

    root_pr_nodes
  }

  if (use_ordering_funct == TRUE) {
    cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds))
  } else {
    cds <- order_cells(cds)
  }


  # Save the pseudotime data and the CDS object to files
  saveRDS(cds, cds_filename)


  # Create and save plots
  gene_plot <- plot_cells(cds = cds, genes = gene_expression, 
                          cell_size = point_size, graph_label_size = 4.5,
                          label_leaves = FALSE, label_branch_points = FALSE,
                          label_cell_groups = FALSE,
                          alpha = 1,
                          label_roots = FALSE)
  gene_plot <- gene_plot + theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_blank(),
    axis.line.x = element_line(size = 1.25, lineend = "round"),
    axis.line.y = element_line(size = 1.25, lineend = "round"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks = element_line(size = 1),
    legend.title = element_text(size = 14, vjust = 1),
    strip.text = element_blank(),
    legend.text = element_text(size = 12)
  ) +
    ggtitle(paste0(gene_title, " Expression"))+
    coord_fixed(ratio = 1, xlim = c(-15, 15), ylim = c(-15, 15))

  ggsave(
    plot = gene_plot, device = plot_type, filename = plot_file,
    width = plot_width, height = plot_height, units = plot_unit, dpi = plot_dpi
  )
  

  pt_graph <- plot_cells(cds,
    color_cells_by = "pseudotime",
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    graph_label_size = 4.5,
    cell_size = point_size,
    alpha = 1
  )
  pt_graph <- pt_graph +
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.line.x = element_line(size = 1.25, lineend = "round"),
          axis.line.y = element_line(size = 1.25, lineend = "round"),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.ticks.length = unit(0.2, "cm"),
          axis.ticks = element_line(size = 1),
          legend.title = element_text(size = 14, vjust = 1),
          legend.text = element_text(size = 12)) +
    ggtitle(paste0(pt_title, " Pseudotime"))+
    scale_color_viridis_c(option = "C") +
    labs(color = "Pseudotime")+
    coord_fixed(ratio = 1, xlim = c(-15, 15), ylim = c(-15, 15))

  ggsave(
    plot = pt_graph, device = plot_type, filename = pt_file,
    width = plot_width, height = plot_height, units = plot_unit, dpi = plot_dpi
  )
  

  # Saving pseudotime data to a dataframe
  pseudotime_data <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  pseudotime_data <- as.data.frame(pseudotime_data)
  colnames(pseudotime_data) <- "Pseudotime"
  pseudotime_contents <- rownames(pseudotime_data)
  pseudotime_data$Samples <- pseudotime_contents

  # Return a list of the relevant objects
  return(list(CDS = cds, Cell_Progression_Graph = gene_plot, PT_Graph = pt_graph, PT_Data = pseudotime_data))
}
