# Name: switchde_calculator.R
# Author: Andrew Willems <awillems@vols.utk.edu>

# Description
# This code calculates the switchde (SDE) metric for gene expression along 
# single-cell RNA seq trajectories and returns a ranking of genes based on their
# switch-like differential expression. The resulting ranking is stored in sde_ranking.

# Input arguments
# - sde: A data frame containing the switchde output.
# - qval: The significance threshold for filtering the switchde output.
# - k: The switchde score of each gene.
# - gene: The gene names corresponding to each switchde score.

# Output
# - sde_ranking: A numeric vector containing the switchde ranking of genes.

# The switchde (SDE) metric
switchde_calculator <- function(denoised_sc = sc,
                                pseudo_time = pt,
                                zero_inflated = FALSE) {
  
  # Loading needed package
  # switchde and tidyverse packages are loaded.
  require(switchde)
  require(tidyverse)
  
  # Actual metric
  # Calculates the Switch-Like Differential Expression (SDE) using the switchde function.
  # Input arguments
  # - denoised_sc: the denoised single-cell RNA-seq data set to analyze.
  # - pseudo_time: the pseudotime information for each cell in the data set.
  # - zero_inflated: a flag indicating if the expression data contains zero-inflated genes.
  sde <- switchde(denoised_sc, as.numeric(pseudo_time$Pseudotime),
                   verbose = TRUE, zero_inflated = zero_inflated)
  
  # Filtering SDE output to just genes that are < 0.05 and then ranking them in
  # decreasing order
  sde_filtered <- filter(sde, qval < 0.05)
  index <- order(abs(sde_filtered$k), decreasing = T)
  sde_rank <- sde_filtered[index,]
  sde_ranking <- sde_rank$k
  names(sde_ranking) <- sde_rank$gene
  sde_ranking<-abs(sde_ranking)/sum(abs(sde_ranking))
  
  # Return object
  return(sde_ranking)
}
