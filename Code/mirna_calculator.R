# Name: mirna_calculator.R
# Author: Andrew Willems <awillems@vols.utk.edu>
# Purpose: Efficiently build miRNA metric outputs

mirna_calculator <- function(ts_org = "Human",
                             ts_version = "8.0",
                             max_mir_targets = 10,
                             cancer_up = TRUE,
                             cancer_type1 = "lung cancer",
                             print_ts_targets = TRUE,
                             mirna_remove = "hsa-miR-129-1-3p",
                             max_mirnas = 807,
                             save_common_mirnas = FALSE,
                             common_mirnas_name = "common_mirnas.csv",
                             save_mirna_genes = TRUE,
                             save_mirna_genes_mat = TRUE,
                             save_mirna_genes_mat_r = TRUE,
                             save_mirna_genes_mat_csv = FALSE,
                             include_sites_score = FALSE,
                             mirna_genes_mat_name = "Outputs/mirna_score_matrix.rds",
                             mirna_ranking_name = "Outputs/mirna_score.rds") {
  # Loading required package
  require(hoardeR)
  require(tidyverse)


  # miRNAs from miRmap
  mirmap_mirnas <- read.csv(
    file = "Data/mirna_data/mirmap_mirnas.csv",
    sep = ","
  )

  if (cancer_up == TRUE) {
    dbdemc_high <- read.csv(
      file = "Data/mirna_data/dbdemc_2_0_high.txt",
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
      file = "Data/mirna_data/dbdemc_2_0_low.txt",
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

  # Finding intersection of the databases to get a set of common miRNAs
  common_mirnas <- common_mirnas[!common_mirnas %in% mirna_remove]
  if (save_common_mirnas == TRUE) {
    write.csv(common_mirnas, file = paste0(
      "Outputs/mirna/common_mirnas/",
      common_mirnas_name
    ))
  }
  if (length(common_mirnas) >= length(common_mirnas[1:max_mirnas])) {
    common_mirnas <- common_mirnas[1:max_mirnas]
  } else {
    print(paste0(
      "There are fewer target miRNAs available than your input.",
      "Using the largest number of common miRNAs for",
      "this submission to TargetScan"
    ))
    common_mirnas <- common_mirnas[1:length(common_mirnas)]
    print(paste0("The number of common miRNAs is: ", length(common_mirnas)))
  }

  if (print_ts_targets == TRUE) {
    print(length(common_mirnas))
  }

  # Submitting all common miRNAs to TargetScan
  mirna_score <- targetScan(
    mirna = common_mirnas,
    species = ts_org,
    release = ts_version,
    maxOut = max_mir_targets
  )
  

  # Collapsing the list of data frames to a single data frame
  mirna_score_df_raw <- bind_rows(mirna_score)

  # Making the scoring matrix that we will update
  mirna_score_mat <- matrix(
    data = 0, nrow = length(unique(mirna_score_df_raw$Ortholog)),
    ncol = length(names(mirna_score)),
    dimnames = list(
      unique(mirna_score_df_raw$Ortholog),
      names(mirna_score)
    )
  )

  mirna_score_df <- as.data.frame(mirna_score_mat)

  # Going through each miRNA-gene pair and updating the score if the
  # respective gene is a target of said miRNA.
  for (m in names(mirna_score)) {
    mirna_score_df[, m] <- ifelse(unique(mirna_score_df_raw$Ortholog) %in% mirna_score[[paste(m)]]$Ortholog, 1, 0)
  }

  if (save_mirna_genes_mat == TRUE & save_mirna_genes_mat_r == TRUE) {
    saveRDS(mirna_score_df, file = mirna_genes_mat_name)
  } else if (save_mirna_genes_mat == TRUE & save_mirna_genes_mat_csv == TRUE) {
    write.csv(mirna_score_df, file = mirna_genes_mat_name)
  }

  # Now calculating the row sums of each gene for total number of
  # the miRNA interactions
  mirna_gene_list <- rowSums(mirna_score_df)
  mirna_ranking <- abs(mirna_gene_list) / sum(abs(mirna_gene_list))
  mirna_ranking <- sort(mirna_ranking, decreasing = TRUE)
  if (save_mirna_genes == TRUE) {
    saveRDS(mirna_ranking, file = mirna_ranking_name)
  }
  # Return object
  return(mirna_ranking)
}
