mirna_calculator <- function(ts_org = "Human",
                             ts_version = "8.0",
                             max_mir_targets = 10,
                             cancer_up = TRUE,
                             status = "up",
                             cancer_type1 = "lung cancer",
                             print_ts_targets = TRUE,
                             mirna_remove = "hsa-miR-129-1-3p",
                             max_mirnas = 808,
                             save_mirna_genes = TRUE,
                             save_mirna_genes_mat = TRUE,
                             include_sites_score = FALSE,
                             save_site_score = FALSE,
                             mirna_ranking_name_rds = "Outputs/lung_mirnas.rds",
                             mirna_genes_mat_name_rds = "Outputs/lung_mirna_score_matrix.rds",
                             mirna_genes_mat_name = "Outputs/lung_mirna_score_matrix.csv",
                             mirna_ranking_name = "Outputs/lung_mirnas.csv",
                             mirna_site_score_name = "Outputs/site_scores.csv") {
  library(furrr)
  library(future)
  library(hoardeR)
  library(retry)

  # Set up the number of cores for parallelization
  plan(multisession, workers = 11)



  # miRNAs from miRmap
  mirmap_mirnas <- read.csv("Data/mirna_data/mirmap_mirnas.csv", sep = ",")

  # Read and filter cancer data
  cancer_data <- read.csv(file = if (cancer_up) "Data/mirna_data/dbdemc_2.0_high.txt" else "Data/mirna_data/dbdemc_2.0_low.txt", sep = "\t")
  if (cancer_up == TRUE && status == "up") {
    cancer_data <- cancer_data %>%
      filter(Cancer.Type == cancer_type1) %>%
      filter(Status == "UP")
  } else if (cancer_up == TRUE && status == "down") {
    cancer_data <- cancer_data %>%
      filter(Cancer.Type == cancer_type1) %>%
      filter(Status == "DOWN")
  }

  # Extract miRNA names
  cancer_mirnas <- unique(cancer_data$miRBase.Update.ID)
  cancer_mirnas <- cancer_mirnas[cancer_mirnas != "unknown"]

  # Common miRNAs
  common_mirnas <- intersect(mirmap_mirnas$mature_name, cancer_mirnas)
  common_mirnas <- common_mirnas[!common_mirnas %in% mirna_remove]

  # Limit number of miRNAs
  if (length(common_mirnas) > max_mirnas) {
    common_mirnas <- common_mirnas[1:max_mirnas]
  } else if (print_ts_targets) {
    cat("There are fewer target miRNAs available than your input. Using the largest number of common miRNAs for this submission to TargetScan\n")
    cat("The number of common miRNAs is:", length(common_mirnas), "\n")
  }

  # TargetScan function with try-catch to help manage no table returned errors
  mirna_targets <- future_map(common_mirnas, function(m) {
    result <- NULL
    tryCatch(
      {
        result <- hoardeR::targetScan(
          mirna = m,
          species = ts_org,
          release = ts_version,
          maxOut = max_mir_targets
        ) %>%
          mutate(mirna_name_final = m)
      },
      error = function(e) {
        message(paste0("Error: ", e$message, ". Moving to next miRNA"))
      }
    )
    result
  })

  # Simplify TargetScan results
  if (include_sites_score) {
    total_list <- lapply(mirna_targets, function(df) {
      if ("Ortholog" %in% names(df)) {
        df <- df %>%
          select(Ortholog, mirna_name_final, consSites, poorlySites) %>%
          rename(gene_list = Ortholog, mirna_list = mirna_name_final) %>%
          mutate(across(c(consSites, poorlySites), ~ as.numeric(gsub("\\*", "", .))), .keep = "unused") %>%
          mutate(across(c(consSites, poorlySites), ~ replace(., is.na(.), 0))) %>%
          mutate(site_score = (consSites / (consSites + poorlySites)) + consSites) %>%
          arrange(desc(site_score))
        df
      } else {
        writeLines("This miRNA is not found in TargetScan. We are continuing on to the next miRNA")
        NULL
      }
    })
  } else {
    total_list <- lapply(mirna_targets, function(df) {
      if ("Ortholog" %in% names(df)) {
        df %>%
          select(Ortholog, mirna_name_final) %>%
          rename(gene_list = Ortholog, mirna_list = mirna_name_final)
      } else {
        writeLines("This miRNA is not found in TargetScan. We are continuing on to the next miRNA")
        NULL
      }
    })
  }



  # Removing any lists in total_list that are NULL (i.e. aren't found in TargetScan)
  total_list <- total_list[!sapply(total_list, is.null)]

  # Function to sum the site scores across unique miRNA-gene pairs
  sum_site_scores <- function(lists) {
    # Create an empty data frame to store the summed site scores
    result <- data.frame(miRNA = character(), gene = character(), site_score = numeric(), stringsAsFactors = FALSE)

    # Loop through each list in the input
    for (lst in lists) {
      # Extract the miRNA, gene, and site score from the list
      miRNA <- lst$mirna_list
      gene <- lst$gene_list
      site_score <- lst$site_score

      # Check if the miRNA-gene pair already exists in the result data frame
      if (any(result$miRNA == miRNA & result$gene == gene)) {
        # If the pair exists, update the site score by summing with the existing value
        result$site_score[result$miRNA == miRNA & result$gene == gene] <- result$site_score[result$miRNA == miRNA & result$gene == gene] + site_score
      } else {
        # If the pair doesn't exist, add a new row to the result data frame
        result <- rbind(result, data.frame(miRNA = miRNA, gene = gene, site_score = site_score, stringsAsFactors = FALSE))
      }
    }

    # Return the resulting data frame
    result <- result[order(result$site_score, decreasing = TRUE), ]
    result
  }

  # Call the function with the input lists
  summed_site_scores <- sum_site_scores(total_list)

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
      if (!is.na(gene)) {
        tryCatch(
          {
            mirna_score[gene, mirna] <- mirna_score[gene, mirna] + 1
          },
          error = function(e) {
            message(paste(
              "Error with gene", gene, "and miRNA", mirna, ":",
              conditionMessage(e)
            ))
          }
        )
      }
    }
  }




  # Save score matrix
  if (save_mirna_genes_mat) {
    write.csv(mirna_score, file = mirna_genes_mat_name)
    saveRDS(mirna_score, file = mirna_genes_mat_name_rds)
  }

  # Calculate and save miRNA ranking
  mirna_ranking <- rowSums(mirna_score) / sum(mirna_score)
  mirna_ranking <- sort(mirna_ranking, decreasing = TRUE)

  if (save_mirna_genes) {
    write.csv(mirna_ranking, file = mirna_ranking_name)
    saveRDS(mirna_ranking, file = mirna_ranking_name_rds)
  }

  if (save_site_score) {
    write.csv(summed_site_scores, file = mirna_site_score_name)
  }

  return(mirna_ranking)
}
