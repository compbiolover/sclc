mirna_calculator <- function(ts_org = "Human",
                             ts_version = "8.0",
                             max_mir_targets = 10,
                             cancer_up = TRUE,
                             cancer_type1 = "lung cancer",
                             print_ts_targets = TRUE,
                             mirna_remove = "hsa-miR-129-1-3p",
                             max_mirnas = 808,
                             save_mirna_genes = TRUE,
                             save_mirna_genes_mat = TRUE,
                             include_sites_score = FALSE,
                             mirna_genes_mat_name = "Outputs/lung_mirna_score_matrix.rds",
                             mirna_ranking_name = "Outputs/lung_mirnas.rds") {
  library(hoardeR)

  # Set up the number of cores for parallelization
  plan(multisession, workers = 11)

  # Timeout error code
  is_timeout <- function(err) {
    inherits(err, "TimeoutError")
  }

  stop_after_n_tries <- function(n) {
    function(err, try_num) {
      if (try_num >= n) {
        stop("Retry limit reached")
      }
      TRUE
    }
  }


  # miRNAs from miRmap
  mirmap_mirnas <- read.csv("Data/mirna_data/mirmap_mirnas.csv", sep = ",")

  # Read and filter cancer data
  cancer_data <- read.csv(file = if (cancer_up) "Data/mirna_data/dbdemc_2.0_high.txt" else "Data/mirna_data/dbdemc_2.0_low.txt", sep = "\t")
  cancer_data <- cancer_data %>%
    filter(Cancer.Type == cancer_type1) %>%
    filter(Status == "UP")

  # Extract miRNA names
  cancer_mirnas <- unique(cancer_data$miRBase.Update.ID)
  cancer_mirnas <- cancer_mirnas[cancer_mirnas != "unknown"]

  # Common miRNAs
  common_mirnas <- intersect(mirmap_mirnas$mature_name, cancer_mirnas)
  common_mirnas <- common_mirnas[!common_mirnas %in% mirna_remove]
  saveRDS(common_mirnas, file = "Outputs/common_mirnas.rds")

  # Limit number of miRNAs
  if (length(common_mirnas) > max_mirnas) {
    common_mirnas <- common_mirnas[1:max_mirnas]
  } else if (print_ts_targets) {
    cat("There are fewer target mirnas available than your input. Using the largest number of common mirnas for this submission to TargetScan\n")
    cat("The number of common mirnas is:", length(common_mirnas), "\n")
  }



  # Get TargetScan results
  mirna_targets <- future_map(common_mirnas, function(m) {
    print(paste0("Processing miRNA: "), m)
    hoardeR::targetScan(
      mirna = m,
      species = ts_org,
      release = ts_version,
      maxOut = max_mir_targets
    ) %>%
      mutate(miRNA_name_final = m)
  })


  # Simplify TargetScan results
  total_list <- lapply(mirna_targets, function(df) {
    if ("Ortholog" %in% names(df)) {
      df %>%
        select(Ortholog, miRNA_name_final) %>%
        rename(gene_list = Ortholog, mirna_list = miRNA_name_final)
    } else {
      writeLines("This miRNA is not found in TargetScan. We are continuing on to the next miRNA")
      NULL
    }
  })
  
  
  # Removing any lists in total_list that are NULL (i.e. aren't found in TargetScan)
  total_list <- total_list[!sapply(total_list, is.null)]

  # Create score matrix
  all_mirs_for_score <- unique(unlist(sapply(total_list, "[[", "mirna_list")))
  all_genes_for_score <- unique(unlist(sapply(total_list, "[[", "gene_list")))
  mirna_score <- matrix(0,
    nrow = length(all_genes_for_score), ncol = length(all_mirs_for_score),
    dimnames = list(all_genes_for_score, all_mirs_for_score)
  )

  # Calculate scores
  for (i in seq_along(total_list)) {
    df <- total_list[[i]]
    for (j in seq_len(nrow(df))) {
      gene <- df[j, "gene_list"]
      mirna <- df[j, "mirna_list"]
      mirna_score[gene, mirna] <- mirna_score[gene, mirna] + 1
    }
  }

  # Save score matrix
  if (save_mirna_genes_mat) {
    saveRDS(mirna_score, file = mirna_genes_mat_name)
  }

  # Calculate and save miRNA ranking
  mirna_ranking <- rowSums(mirna_score) / sum(mirna_score)
  mirna_ranking <- sort(mirna_ranking, decreasing = TRUE)

  if (save_mirna_genes) {
    saveRDS(mirna_ranking, file = mirna_ranking_name)
  }

  return(mirna_ranking)
}
