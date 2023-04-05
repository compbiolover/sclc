# Name: mirna_calculator.R
# Author: Andrew Willems <awillems@vols.utk.edu>
# Purpose: Efficiently build miRNA metric outputs

# mirna_calculator
mirna_calculator <- function(ts_org = "Human",
                             ts_version = "8.0",
                             max_mir_targets = 10,
                             cancer_up = TRUE,
                             cancer_type1 = "lung cancer",
                             print_ts_targets = TRUE,
                             mirna_remove = "hsa-miR-129-1-3p",
                             max_mirnas = 807,
                             save_mirna_genes = TRUE,
                             save_mirna_genes_mat = TRUE,
                             include_sites_score = FALSE,
                             mirna_genes_mat_name = "Outputs/lung_mirna_score_matrix.rds",
                             mirna_ranking_name = "Outputs/lung_mirnas.rds") {

  library(hoardr)
  # miRNAs from miRmap
  mirmap_mirnas <- read.csv(
    file = "Data/miRNA-data/mirmap_mirnas.csv",
    sep = ","
  )

  if (cancer_up == TRUE) {
    dbdemc_high <- read.csv(
      file = "Data/miRNA-data/dbDEMC-2.0-high.txt",
      sep = "\t"
    )

    # Filtering to just the miRNAs associated with a particular type of cancer
    dbdemc_high <- filter(
      dbdemc_high,
      Cancer.Type == cancer_type1
    )
    dbdemc_high_mirnas <- subset(dbdemc_high, select = miRBase.Update.ID)
    if ("unknown" %in% dbdemc_high_mirnas$miRBase.Update.ID == TRUE) {
      dbdemc_high_mirnas <- filter(
        dbdemc_high_mirnas,
        miRBase.Update.ID != "unknown"
      )
      dbdemc_high_mirnas <- as.vector(dbdemc_high_mirnas)
    } else {
      dbdemc_high_mirnas <- as.vector(dbdemc_high_mirnas)
    }

    # Common miRNAs
    common_mirnas <- intersect(
      mirmap_mirnas$mature_name,
      dbdemc_high_mirnas$miRBase.Update.ID
    )
  } else {
    dbdemc_low <- read.csv(
      file = "Data/miRNA-data/dbDEMC-2.0-low.txt",
      sep = "\t"
    )
    colnames(dbdemc_low)[1] <- "miRNA.ID"

    # Filtering to just the miRNAs associated with a particular type of cancer
    dbdemc_low <- filter(
      dbdemc_low,
      Cancer_Type == cancer_type1
    )
    dbdemc_low_mirnas <- subset(dbdemc_low, select = miRBase_update)
    dbdemc_low_mirnas <- as.vector(dbdemc_low_mirnas)

    # Common miRNAs between databases
    common_mirnas <- intersect(
      mirmap_mirnas$mature_name,
      dbdemc_low_mirnas$miRBase_update
    )
  }

  # Now submitting these miRNAs to TargetScan to get genes to make a gene list
  # for the third metric testing to see if all of the miRNAs exist in TargetScan
  # before getting all submitted at once
  common_mirnas <- common_mirnas[!common_mirnas %in% mirna_remove]
  saveRDS(common_mirnas, file = "Outputs/common_mirnas.rds")
  if (length(common_mirnas) >= length(common_mirnas[1:max_mirnas])) {
    common_mirnas <- common_mirnas[1:max_mirnas]
  } else {
    print("There are fewer target mirnas available than your input. Using the largest number of common mirnas for this submission to TargetScan")
    common_mirnas <- common_mirnas[1:length(common_mirnas)]
    print(paste0("The number of common mirnas is: ", length(common_mirnas)))
  }

  if (print_ts_targets == TRUE) {
    print(length(common_mirnas))
  }
  my_num <- 1
  mirna_targets <- vector(mode = "list", length = length(common_mirnas))
  for (m in common_mirnas[1:length(common_mirnas)]) {
    percent_done <- (my_num / length(common_mirnas)) * 100
    print(paste0("Currently ", round(percent_done, digits = 2), "% done"))
    print(paste0("We are currently on miRNA: ", m))
    current_target <- hoardeR::targetScan(
      mirna = common_mirnas[my_num],
      species = ts_org,
      release = ts_version,
      maxOut = max_mir_targets
    )

    mirna_name <- m
    mirna_name_final <- rep(mirna_name, times = length(current_target$Ortholog))
    current_target <- cbind(current_target, mirna_name_final)
    mirna_targets[[m]] <- current_target
    my_num <- my_num + 1
  }



  # Simplifying the output of the TargetScan commands to
  # just the Gene name and the miRNA columns
  counter <- 1
  total_list <- list()
  for (i in mirna_targets) {
    current_df <- i
    gene_list <- current_df$Ortholog
    mirna_list <- current_df$miRNA_name_final
    simple_list <- cbind(gene_list, mirna_list)
    total_list[[counter]] <- simple_list
    counter <- counter + 1
  }

  counter <- 1
  all_mirs_for_score <- list()
  all_genes_for_score <- list()
  for (l in total_list) {
    current_df <- l
    current_mirs <- current_df[, "mirna_list"]
    all_mirs_for_score[[counter]] <- unique(current_mirs)
    current_genes <- current_df[, "gene_list"]
    all_genes_for_score[[counter]] <- current_genes
    counter <- counter + 1
  }

  all_mirs_for_score <- unlist(all_mirs_for_score)
  all_genes_for_score <- unlist(all_genes_for_score)
  all_genes_for_score_unique <- unique(all_genes_for_score)

  mirna_score <- matrix(
    data = 0, nrow = length(all_genes_for_score_unique),
    ncol = length(all_mirs_for_score),
    dimnames = list(
      all_genes_for_score_unique,
      all_mirs_for_score
    )
  )
  mirna_score <- as.data.frame(mirna_score)

  # Checking to see for each miRNA (column name) if it interacts with a particular row.
  # If it does it get a plus one to that cell. If it does not it moves to next cell
  for (i in rownames(mirna_score)) {
    for (x in mirna_targets) {
      current_df <- x
      if (i %in% current_df$Ortholog) {
        mirna_to_add <- unique(current_df$miRNA_name_final)
        mirna_score[i, mirna_to_add] <- mirna_score[i, mirna_to_add] + 1
      }
    }
  }

  if (save_mirna_genes_mat == TRUE) {
    saveRDS(mirna_score, file = mirna_genes_mat_name)
  }

  # Now calculating the row sums of each gene for total number of
  # the miRNA interactions
  mirna_gene_list <- rowSums(mirna_score)
  mirna_ranking <- abs(mirna_gene_list) / sum(abs(mirna_gene_list))
  mirna_ranking <- sort(mirna_ranking, decreasing = TRUE)
  if (save_mirna_genes == TRUE) {
    saveRDS(mirna_ranking, file = mirna_ranking_name)
  }
  # Return object
  return(mirna_ranking)
}
