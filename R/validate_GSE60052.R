#!/usr/bin/env Rscript
# ============================================================================
# External Validation: 4-Gene miRMS Signature on GSE60052
# ============================================================================
#
# Run this script locally where you have network access:
#   Rscript R/validate_GSE60052.R
#
# Or source it from within main.Rmd after placing the downloaded data files
# in data/external/GSE60052/
#
# Requirements: GEOquery, Biobase, survival, survminer, timeROC, ggplot2
# ============================================================================

# ── Install missing packages ────────────────────────────────────────────────
required_pkgs <- c("GEOquery", "Biobase", "survival", "survminer",
                   "timeROC", "ggplot2", "dplyr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("GEOquery", "Biobase")) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}

library(GEOquery)
library(Biobase)
library(survival)
library(survminer)
library(timeROC)
library(ggplot2)
library(dplyr)

# ── Configuration ───────────────────────────────────────────────────────────
SIGNATURE_GENES <- c("TFRC", "FAM83F", "DLK1", "GNG13")
SIGNATURE_COEFS <- c(-0.9526172, -0.0854650, -0.0006146, -0.2557644)
names(SIGNATURE_COEFS) <- SIGNATURE_GENES

OUTPUT_DIR <- file.path(getwd(), "Figures", "panels")
DATA_CACHE <- file.path(getwd(), "data", "external", "GSE60052")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_CACHE, recursive = TRUE, showWarnings = FALSE)

cat("=" , rep("=", 69), "\n", sep = "")
cat("External Validation: 4-Gene miRMS Signature on GSE60052\n")
cat("Signature: TFRC, FAM83F, DLK1, GNG13\n")
cat("=" , rep("=", 69), "\n\n", sep = "")

# ── Step 1: Download GSE60052 ───────────────────────────────────────────────
cat("STEP 1: Downloading GSE60052 from GEO...\n")

# GSE60052 is RNA-seq (GPL11154 = Illumina HiSeq 2000)
# GEOquery may return the series matrix (processed expression) and phenoData
gse <- getGEO("GSE60052", destdir = DATA_CACHE, GSEMatrix = TRUE)

# If multiple platforms, select the right one
if (length(gse) > 1) {
  cat("  Multiple platforms found:", names(gse), "\n")
  cat("  Selecting first ExpressionSet...\n")
}
eset <- gse[[1]]

cat("  Samples:", ncol(eset), "\n")
cat("  Features:", nrow(eset), "\n")
cat("  Platform:", annotation(eset), "\n\n")

# ── Step 2: Extract phenotype/clinical data ─────────────────────────────────
cat("STEP 2: Extracting clinical data...\n")

pheno <- pData(eset)
cat("  Phenotype columns:\n")
print(names(pheno))
cat("\n")

# Print all characteristics to find survival columns
char_cols <- grep("characteristics", names(pheno), value = TRUE, ignore.case = TRUE)
cat("  Characteristics columns:\n")
for (cc in char_cols) {
  cat("    ", cc, ":", head(unique(pheno[[cc]]), 5), "\n")
}
cat("\n")

# Parse characteristics_ch1 fields into structured columns
# GSE60052 typically has: tissue, age, gender, smoking, stage, survival_time,
#                         vital_status in characteristics_ch1.*
parse_characteristic <- function(x) {
  if (is.na(x) || x == "") return(NA)
  parts <- strsplit(as.character(x), ":\\s*")[[1]]
  if (length(parts) >= 2) return(parts[2])
  return(x)
}

# Extract key and value from each characteristic column
clinical_df <- data.frame(
  sample_id = rownames(pheno),
  stringsAsFactors = FALSE
)

for (cc in char_cols) {
  # Get the key name from the first non-NA entry
  first_val <- na.omit(pheno[[cc]])[1]
  key <- trimws(strsplit(as.character(first_val), ":")[[1]][1])
  key <- tolower(gsub("[^a-zA-Z0-9]", "_", key))

  clinical_df[[key]] <- sapply(pheno[[cc]], parse_characteristic)
}

cat("  Parsed clinical columns:", names(clinical_df), "\n\n")

# ── Step 3: Identify survival columns ───────────────────────────────────────
cat("STEP 3: Identifying survival data...\n")

# Look for time/survival/status columns
surv_cols <- grep("surv|time|status|vital|alive|dead|os|follow|event|month|day",
                  names(clinical_df), value = TRUE, ignore.case = TRUE)
cat("  Potential survival columns:", surv_cols, "\n")

# Print unique values to understand encoding
for (sc in surv_cols) {
  cat("    ", sc, ": ", paste(head(unique(clinical_df[[sc]]), 10), collapse = ", "), "\n")
}
cat("\n")

# If no survival columns found in characteristics, check supplementary files
if (length(surv_cols) == 0) {
  cat("  WARNING: No survival columns found in GEO metadata.\n")
  cat("  Checking supplementary files from the original paper...\n")
  cat("  Paper: Jiang et al., PLOS Genetics 2016 (10.1371/journal.pgen.1005895)\n")
  cat("  You may need to download S1 Table from the paper manually.\n\n")

  # Try to download supplementary files
  supp_files <- tryCatch(
    getGEOSuppFiles("GSE60052", baseDir = DATA_CACHE),
    error = function(e) {
      cat("  Could not auto-download supplementary files:", e$message, "\n")
      NULL
    }
  )

  if (!is.null(supp_files)) {
    cat("  Supplementary files downloaded:\n")
    print(rownames(supp_files))
  }
}

# ── Step 4: Extract expression data for signature genes ─────────────────────
cat("STEP 4: Extracting expression data for signature genes...\n")

expr_mat <- exprs(eset)
cat("  Expression matrix dimensions:", dim(expr_mat), "\n")
cat("  Expression value range:", range(expr_mat, na.rm = TRUE), "\n")
cat("  First few row IDs:", head(rownames(expr_mat)), "\n\n")

# Check if row IDs are gene symbols or probe IDs
# For RNA-seq data (GPL11154), rows may already be gene symbols
# For microarray, we may need to map probe IDs to gene symbols

# Try direct gene symbol matching first
gene_matches <- list()
for (gene in SIGNATURE_GENES) {
  matches <- grep(paste0("^", gene, "$"), rownames(expr_mat),
                  ignore.case = TRUE, value = TRUE)
  if (length(matches) > 0) {
    gene_matches[[gene]] <- matches[1]
    cat("  ", gene, ": FOUND as", matches[1], "\n")
  } else {
    # Try partial match
    partial <- grep(gene, rownames(expr_mat), ignore.case = TRUE, value = TRUE)
    if (length(partial) > 0) {
      gene_matches[[gene]] <- partial[1]
      cat("  ", gene, ": partial match ->", partial[1], "\n")
    } else {
      cat("  ", gene, ": NOT FOUND in row names\n")
    }
  }
}

# If genes not found by symbol, try feature annotation table
if (length(gene_matches) < length(SIGNATURE_GENES)) {
  cat("\n  Attempting probe-to-gene mapping via feature data...\n")
  feat <- fData(eset)
  cat("  Feature data columns:", names(feat), "\n")

  # Look for a gene symbol column
  symbol_cols <- grep("symbol|gene.*name|gene_assignment|gene_id|genename",
                      names(feat), value = TRUE, ignore.case = TRUE)
  cat("  Gene symbol columns:", symbol_cols, "\n")

  if (length(symbol_cols) > 0) {
    sym_col <- symbol_cols[1]
    for (gene in SIGNATURE_GENES) {
      if (is.null(gene_matches[[gene]])) {
        probe_idx <- grep(paste0("\\b", gene, "\\b"), feat[[sym_col]],
                          ignore.case = TRUE)
        if (length(probe_idx) > 0) {
          # If multiple probes, take the one with highest variance
          if (length(probe_idx) > 1) {
            vars <- apply(expr_mat[probe_idx, , drop = FALSE], 1, var, na.rm = TRUE)
            probe_idx <- probe_idx[which.max(vars)]
          }
          gene_matches[[gene]] <- rownames(expr_mat)[probe_idx]
          cat("  ", gene, ": mapped via", sym_col, "->", rownames(expr_mat)[probe_idx], "\n")
        }
      }
    }
  }
}

cat("\n  Final gene mapping:\n")
for (gene in SIGNATURE_GENES) {
  if (!is.null(gene_matches[[gene]])) {
    cat("    ", gene, "->", gene_matches[[gene]], "\n")
  } else {
    cat("    ", gene, "-> MISSING\n")
  }
}

found_genes <- names(gene_matches)
n_found <- length(found_genes)

if (n_found < 3) {
  stop("Fewer than 3 signature genes found. Cannot proceed with validation.")
}
if (n_found < 4) {
  cat("\n  WARNING: Only", n_found, "of 4 genes found.",
      "Proceeding with partial signature.\n")
}

# ── Step 5: Build expression matrix for signature genes ─────────────────────
cat("\nSTEP 5: Building expression data for validation...\n")

# Extract expression for found genes
expr_sig <- t(expr_mat[unlist(gene_matches), , drop = FALSE])
colnames(expr_sig) <- names(gene_matches)

cat("  Signature expression matrix:", nrow(expr_sig), "samples x",
    ncol(expr_sig), "genes\n")

# Check if data needs log transformation
max_val <- max(expr_sig, na.rm = TRUE)
if (max_val > 30) {
  cat("  Max expression value:", max_val, "— appears to be raw counts/FPKM.\n")
  cat("  Applying log2(x + 1) transformation.\n")
  expr_sig <- log2(expr_sig + 1)
  log_transformed <- TRUE
} else {
  cat("  Max expression value:", max_val, "— appears already log-transformed.\n")
  log_transformed <- FALSE
}

cat("  Expression summary after transformation:\n")
print(summary(expr_sig))

# ── Step 6: Match clinical and expression data ──────────────────────────────
cat("\nSTEP 6: Matching clinical and expression data...\n")

# Build the survival + expression data frame
validation_df <- data.frame(
  sample_id = rownames(expr_sig),
  expr_sig,
  stringsAsFactors = FALSE
)

# Merge with clinical data
validation_df <- merge(validation_df, clinical_df, by = "sample_id")

cat("  Merged samples:", nrow(validation_df), "\n")
cat("  Columns:", names(validation_df), "\n\n")

# ── Step 7: Identify and prepare survival variables ─────────────────────────
cat("STEP 7: Preparing survival variables...\n")

# Attempt to automatically detect survival time and status columns
# Common patterns in GEO SCLC datasets:
time_candidates <- grep("surv.*time|os.*time|overall.*surv|follow.*up|time.*month|time.*day|survival_time|os_months",
                        names(validation_df), value = TRUE, ignore.case = TRUE)
status_candidates <- grep("status|vital|alive|dead|event|censor|os_status",
                          names(validation_df), value = TRUE, ignore.case = TRUE)

cat("  Time column candidates:", time_candidates, "\n")
cat("  Status column candidates:", status_candidates, "\n")

# If we found survival columns, proceed
if (length(time_candidates) > 0 && length(status_candidates) > 0) {
  time_col <- time_candidates[1]
  status_col <- status_candidates[1]

  cat("  Using time:", time_col, "\n")
  cat("  Using status:", status_col, "\n")

  # Convert to numeric
  validation_df$os_time <- as.numeric(validation_df[[time_col]])
  status_raw <- validation_df[[status_col]]

  # Handle various status encodings
  if (all(status_raw %in% c("0", "1", 0, 1), na.rm = TRUE)) {
    validation_df$os_status <- as.numeric(status_raw)
  } else if (any(grepl("dead|deceased|1", status_raw, ignore.case = TRUE))) {
    validation_df$os_status <- ifelse(
      grepl("dead|deceased|1", status_raw, ignore.case = TRUE), 1, 0)
  } else if (any(grepl("alive|living|0", status_raw, ignore.case = TRUE))) {
    validation_df$os_status <- ifelse(
      grepl("alive|living|0", status_raw, ignore.case = TRUE), 0, 1)
  } else {
    cat("  WARNING: Cannot parse status encoding:", unique(status_raw), "\n")
    cat("  Please set validation_df$os_status manually.\n")
  }
} else {
  cat("\n  === MANUAL INTERVENTION NEEDED ===\n")
  cat("  Could not auto-detect survival columns.\n")
  cat("  Available columns:", names(validation_df), "\n")
  cat("  Unique values per column:\n")
  for (col in names(clinical_df)[-1]) {
    cat("    ", col, ":", paste(head(unique(validation_df[[col]]), 5), collapse = ", "), "\n")
  }
  cat("\n  Please set os_time and os_status columns manually and re-run from Step 8.\n")
  cat("  Saving intermediate data to:", file.path(DATA_CACHE, "validation_df_intermediate.rds"), "\n")
  saveRDS(validation_df, file.path(DATA_CACHE, "validation_df_intermediate.rds"))

  # Also save a CSV for easy inspection
  write.csv(validation_df, file.path(DATA_CACHE, "validation_df_intermediate.csv"),
            row.names = FALSE)
  cat("  Also saved as CSV for inspection.\n")
  cat("  After setting os_time and os_status, source R/validate_GSE60052_step2.R\n")
}

# ── Step 8: Compute risk scores ─────────────────────────────────────────────
cat("\nSTEP 8: Computing risk scores...\n")

# Filter to patients with complete survival data
surv_complete <- validation_df %>%
  filter(!is.na(os_time) & !is.na(os_status) & os_time > 0)

cat("  Patients with survival data:", nrow(surv_complete), "\n")

if (nrow(surv_complete) < 10) {
  cat("  WARNING: Very few patients with survival data. Results may not be reliable.\n")
}

# Compute risk scores using the trained coefficients
used_coefs <- SIGNATURE_COEFS[found_genes]
surv_complete$risk_score <- rowSums(
  sapply(found_genes, function(g) surv_complete[[g]] * used_coefs[g])
)

# Stratify by median
median_risk <- median(surv_complete$risk_score)
surv_complete$risk <- ifelse(surv_complete$risk_score > median_risk, "high", "low")
surv_complete$risk <- factor(surv_complete$risk, levels = c("low", "high"))

cat("  Median risk score:", round(median_risk, 4), "\n")
cat("  High risk:", sum(surv_complete$risk == "high"), "\n")
cat("  Low risk:", sum(surv_complete$risk == "low"), "\n")
cat("  Events (deaths):", sum(surv_complete$os_status == 1), "\n")

# ── Step 9: Kaplan-Meier analysis ───────────────────────────────────────────
cat("\nSTEP 9: Kaplan-Meier survival analysis...\n")

km_fit <- survfit(Surv(os_time, os_status) ~ risk, data = surv_complete)
logrank <- survdiff(Surv(os_time, os_status) ~ risk, data = surv_complete)
lr_pval <- 1 - pchisq(logrank$chisq, df = 1)

cat("  Log-rank p-value:", format.pval(lr_pval, digits = 4), "\n")
cat("  Median survival (high risk):",
    summary(km_fit)$table["risk=high", "median"], "\n")
cat("  Median survival (low risk):",
    summary(km_fit)$table["risk=low", "median"], "\n")

# KM Plot
km_plot <- ggsurvplot(
  fit = km_fit,
  data = surv_complete,
  conf.int = TRUE,
  pval = TRUE,
  pval.size = 5,
  risk.table = TRUE,
  palette = c("blue", "red"),
  censor.shape = "|",
  surv.median.line = "hv",
  size = 1.5,
  censor.size = 4.5,
  xlab = "Time (days)",
  ylab = "Overall Survival Probability",
  legend.title = "Risk Group",
  title = "External Validation: GSE60052 (Chinese SCLC, n=48)",
  legend.labs = c("Low risk", "High risk"),
  ggtheme = theme(
    axis.line = element_line(color = "black", linewidth = 1.25, lineend = "round"),
    axis.title = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.ticks = element_line(linewidth = 1.25),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    panel.background = element_blank()
  )
)

ggsave(km_plot$plot,
       filename = "external_validation_GSE60052_KM.png",
       path = OUTPUT_DIR, dpi = 600, width = 8, height = 8, units = "in")
ggsave(km_plot$plot,
       filename = "external_validation_GSE60052_KM.svg",
       path = OUTPUT_DIR, dpi = 600, width = 8, height = 8, units = "in")

cat("  KM plot saved to:", OUTPUT_DIR, "\n")

# ── Step 10: Cox PH model & C-index ────────────────────────────────────────
cat("\nSTEP 10: Cox PH model and C-index...\n")

cox_model <- coxph(Surv(os_time, os_status) ~ risk_score, data = surv_complete)
c_index <- summary(cox_model)$concordance[1]
c_index_se <- summary(cox_model)$concordance[2]

cat("  C-index:", round(c_index, 4), "\n")
cat("  C-index SE:", round(c_index_se, 4), "\n")
cat("  C-index 95% CI: [",
    round(c_index - 1.96 * c_index_se, 4), ",",
    round(c_index + 1.96 * c_index_se, 4), "]\n")

cat("\n  Cox model summary:\n")
print(summary(cox_model))

# ── Step 11: Time-dependent ROC ─────────────────────────────────────────────
cat("\nSTEP 11: Time-dependent ROC analysis...\n")

# Check max follow-up to set appropriate timepoints
max_time <- max(surv_complete$os_time, na.rm = TRUE)
cat("  Max follow-up time:", max_time, "\n")

# Set timepoints based on available follow-up
# GSE60052 survival is typically in months — check scale
if (max_time < 100) {
  cat("  Time appears to be in months.\n")
  roc_times <- c(12, 24, 36)  # 1, 2, 3 years in months
  roc_labels <- c("1-year", "2-year", "3-year")
  time_unit <- "months"
} else {
  cat("  Time appears to be in days.\n")
  roc_times <- c(365, 730, 1095)  # 1, 2, 3 years in days
  roc_labels <- c("1-year", "2-year", "3-year")
  time_unit <- "days"
}

# Only use timepoints within follow-up range
roc_times <- roc_times[roc_times < max_time * 0.9]

if (length(roc_times) > 0) {
  troc <- tryCatch({
    timeROC(
      T = surv_complete$os_time,
      delta = surv_complete$os_status,
      marker = -surv_complete$risk_score,  # negative because lower score = higher risk
      cause = 1,
      times = roc_times,
      iid = TRUE
    )
  }, error = function(e) {
    cat("  Time-dependent ROC failed:", e$message, "\n")
    NULL
  })

  if (!is.null(troc)) {
    for (i in seq_along(roc_times)) {
      ci <- confint(troc, parm = NULL)
      cat("  ", roc_labels[i], "AUC:", round(troc$AUC[i + 1], 4), "\n")
    }
  }
} else {
  cat("  Follow-up too short for time-dependent ROC at standard timepoints.\n")
}

# ── Step 12: Summary ───────────────────────────────────────────────────────
cat("\n")
cat("=" , rep("=", 69), "\n", sep = "")
cat("VALIDATION SUMMARY\n")
cat("=" , rep("=", 69), "\n", sep = "")
cat("Cohort:          GSE60052 (Jiang et al., PLOS Genetics 2016)\n")
cat("Population:      Chinese SCLC\n")
cat("Patients:       ", nrow(surv_complete), "\n")
cat("Events:         ", sum(surv_complete$os_status == 1), "\n")
cat("Genes found:    ", paste(found_genes, collapse = ", "), "\n")
cat("Log-transform:  ", ifelse(log_transformed, "Applied (log2(x+1))", "Not needed"), "\n")
cat("─────────────────────────────────────────────\n")
cat("C-index:        ", round(c_index, 4), " (SE:", round(c_index_se, 4), ")\n")
cat("Log-rank p:     ", format.pval(lr_pval, digits = 4), "\n")
if (!is.null(troc)) {
  for (i in seq_along(roc_times)) {
    cat(roc_labels[i], "AUC:    ", round(troc$AUC[i + 1], 4), "\n")
  }
}
cat("─────────────────────────────────────────────\n")

# Interpretation
if (c_index > 0.65 && lr_pval < 0.05) {
  cat("\nINTERPRETATION: The 4-gene signature shows significant risk\n")
  cat("stratification in the independent GSE60052 cohort.\n")
  cat("C-index >0.65 with p<0.05 supports external validity.\n")
} else if (c_index > 0.60) {
  cat("\nINTERPRETATION: Modest discriminative ability (C-index 0.60-0.65).\n")
  cat("The signature shows some prognostic value but may need refinement.\n")
} else {
  cat("\nINTERPRETATION: Poor discriminative ability (C-index <0.60).\n")
  cat("The signature does not generalize well to this cohort.\n")
  cat("Consider: different expression platforms, population effects,\n")
  cat("or potential overfitting in the training cohort.\n")
}

# Save all results
results <- list(
  cohort = "GSE60052",
  n_patients = nrow(surv_complete),
  n_events = sum(surv_complete$os_status == 1),
  genes_found = found_genes,
  c_index = c_index,
  c_index_se = c_index_se,
  logrank_p = lr_pval,
  cox_model = cox_model,
  km_fit = km_fit,
  time_roc = if (exists("troc")) troc else NULL,
  surv_data = surv_complete,
  log_transformed = log_transformed
)

saveRDS(results, file.path(DATA_CACHE, "validation_results_GSE60052.rds"))
cat("\nResults saved to:", file.path(DATA_CACHE, "validation_results_GSE60052.rds"), "\n")
