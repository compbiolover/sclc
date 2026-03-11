# Name: download_mirna_data.R
# Author: Andrew Willems <awillems@vols.utk.edu>
# Purpose: Download and organize miRNA database files required by mirna_calculator.R.
#
# Usage:
#   source("R/download_mirna_data.R")
#   download_mirna_data()            # downloads everything
#   download_mirna_data("targetscan") # downloads only TargetScan files
#
# After running, your Data/mirna_data/ directory should contain:
#   From TargetScan 8.0 (https://www.targetscan.org/vert_80/):
#     - Summary_Counts.default_predictions.txt
#     - miR_Family_Info.txt
#     - Predicted_Targets_Context_Scores.default_predictions.txt (optional)
#
#   From miRDB 6.0 (https://mirdb.org/download.html):
#     - miRDB_v6.0_prediction_result.txt
#
#   From miRTarBase (https://mirtarbase.cuhk.edu.cn/):
#     - miRTarBase.csv
#
#   From miRmap (https://mirmap.ezlab.org/):
#     - mirmap_mirnas.csv
#
#   From dbDEMC 3.0 (https://www.biosino.org/dbDEMC/):
#     - dbdemc_3.0_high.txt  (upregulated miRNAs)
#     - dbdemc_3.0_low.txt   (downregulated miRNAs)
#
# Note: GEO data (GSE19945) is downloaded automatically at runtime by
#       mirna_calculator() via GEOquery — no manual setup needed for that.


#' Download and set up miRNA database files
#'
#' @param which Character vector specifying which databases to download.
#'   Options: "targetscan", "mirdb", "mirtarbase", "mirmap", "dbdemc", or "all".
#' @param data_dir Directory to save files into. Defaults to "Data/mirna_data".
#' @param overwrite If TRUE, re-download files that already exist.
#' @param verbose Print progress messages.
#'
#' @details
#' Some databases require manual download due to license agreements or lack of
#' direct download URLs. This function will download what it can automatically
#' and print instructions for the rest.
download_mirna_data <- function(which = "all",
                                data_dir = "Data/mirna_data",
                                overwrite = FALSE,
                                verbose = TRUE) {
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
    if (verbose) message("Created directory: ", data_dir)
  }

  targets <- if ("all" %in% which) {
    c("targetscan", "mirdb", "mirtarbase", "mirmap", "dbdemc")
  } else {
    which
  }

  # ---------------------------------------------------------------------------
  # TargetScan 8.0
  # ---------------------------------------------------------------------------
  if ("targetscan" %in% targets) {
    if (verbose) message("\n=== TargetScan 8.0 ===")

    ts_files <- list(
      list(
        url = "https://www.targetscan.org/vert_80/vert_80_data_download/Summary_Counts.default_predictions.txt.zip",
        dest = "Summary_Counts.default_predictions.txt",
        desc = "Default CNN-based predictions (required)"
      ),
      list(
        url = "https://www.targetscan.org/vert_80/vert_80_data_download/miR_Family_Info.txt.zip",
        dest = "miR_Family_Info.txt",
        desc = "miRNA family-to-name mapping (required)"
      ),
      list(
        url = "https://www.targetscan.org/vert_80/vert_80_data_download/Predicted_Targets_Context_Scores.default_predictions.txt.zip",
        dest = "Predicted_Targets_Context_Scores.default_predictions.txt",
        desc = "Per-site context++ scores (optional, improves weighting)"
      )
    )

    for (f in ts_files) {
      dest_path <- file.path(data_dir, f$dest)
      if (file.exists(dest_path) && !overwrite) {
        if (verbose) message("  Already exists: ", f$dest)
        next
      }

      if (verbose) message("  Downloading: ", f$desc)
      zip_path <- tempfile(fileext = ".zip")
      tryCatch({
        download.file(f$url, zip_path, mode = "wb", quiet = !verbose)
        unzip(zip_path, exdir = data_dir)
        if (verbose) message("  -> Saved: ", dest_path)
      }, error = function(e) {
        message("  ! Download failed: ", e$message)
        message("  Manual download: https://www.targetscan.org/cgi-bin/targetscan/data_download.vert80.cgi")
      }, finally = {
        unlink(zip_path)
      })
    }
  }

  # ---------------------------------------------------------------------------
  # miRDB 6.0
  # ---------------------------------------------------------------------------
  if ("mirdb" %in% targets) {
    if (verbose) message("\n=== miRDB 6.0 ===")

    dest_path <- file.path(data_dir, "miRDB_v6.0_prediction_result.txt")
    if (file.exists(dest_path) && !overwrite) {
      if (verbose) message("  Already exists: miRDB_v6.0_prediction_result.txt")
    } else {
      mirdb_url <- "https://mirdb.org/download/miRDB_v6.0_prediction_result.txt.gz"
      if (verbose) message("  Downloading miRDB predictions...")
      gz_path <- tempfile(fileext = ".gz")
      tryCatch({
        download.file(mirdb_url, gz_path, mode = "wb", quiet = !verbose)
        con <- gzfile(gz_path, "rb")
        raw <- readLines(con)
        close(con)
        writeLines(raw, dest_path)
        if (verbose) message("  -> Saved: ", dest_path)
      }, error = function(e) {
        message("  ! Download failed: ", e$message)
        message("  Manual download: https://mirdb.org/download.html")
        message("  Save as: ", dest_path)
      }, finally = {
        unlink(gz_path)
      })
    }
  }

  # ---------------------------------------------------------------------------
  # miRTarBase — requires manual download
  # ---------------------------------------------------------------------------
  if ("mirtarbase" %in% targets) {
    if (verbose) message("\n=== miRTarBase ===")

    dest_path <- file.path(data_dir, "miRTarBase.csv")
    if (file.exists(dest_path) && !overwrite) {
      if (verbose) message("  Already exists: miRTarBase.csv")
    } else {
      message("  miRTarBase requires manual download:")
      message("  1. Go to: https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/php/download.php")
      message("  2. Download the 'miRTarBase' Excel/CSV file for Homo sapiens")
      message("  3. Save as: ", dest_path)
    }
  }

  # ---------------------------------------------------------------------------
  # miRmap — requires manual download
  # ---------------------------------------------------------------------------
  if ("mirmap" %in% targets) {
    if (verbose) message("\n=== miRmap ===")

    dest_path <- file.path(data_dir, "mirmap_mirnas.csv")
    if (file.exists(dest_path) && !overwrite) {
      if (verbose) message("  Already exists: mirmap_mirnas.csv")
    } else {
      message("  miRmap requires manual download:")
      message("  1. Go to: https://mirmap.ezlab.org/")
      message("  2. Export the miRNA list as CSV")
      message("  3. Save as: ", dest_path)
    }
  }

  # ---------------------------------------------------------------------------
  # dbDEMC 3.0 — requires manual download
  # ---------------------------------------------------------------------------
  if ("dbdemc" %in% targets) {
    if (verbose) message("\n=== dbDEMC 3.0 (only needed if use_geo = FALSE) ===")

    high_path <- file.path(data_dir, "dbdemc_3.0_high.txt")
    low_path <- file.path(data_dir, "dbdemc_3.0_low.txt")

    if (file.exists(high_path) && file.exists(low_path) && !overwrite) {
      if (verbose) message("  Already exists: dbdemc files")
    } else {
      message("  dbDEMC 3.0 requires manual download:")
      message("  1. Go to: https://www.biosino.org/dbDEMC/")
      message("  2. Search for 'lung cancer'")
      message("  3. Download upregulated miRNAs -> save as: ", high_path)
      message("  4. Download downregulated miRNAs -> save as: ", low_path)
      message("  Note: dbDEMC 3.0 has 2x more data entries than 2.0 (3,268 unique DEMs")
      message("  across 40 cancer types) and may include fold change/Q-value columns.")
      message("  These files are only needed for the legacy approach (use_geo = FALSE).")
      message("  The GEO-based approach (default) does not require these files.")
    }
  }

  # ---------------------------------------------------------------------------
  # Summary
  # ---------------------------------------------------------------------------
  if (verbose) {
    message("\n=== Setup Summary ===")
    expected <- c(
      "Summary_Counts.default_predictions.txt",
      "miR_Family_Info.txt",
      "miRDB_v6.0_prediction_result.txt",
      "miRTarBase.csv",
      "mirmap_mirnas.csv"
    )
    optional <- c(
      "Predicted_Targets_Context_Scores.default_predictions.txt",
      "dbdemc_3.0_high.txt",
      "dbdemc_3.0_low.txt"
    )

    message("Required files:")
    for (f in expected) {
      status <- if (file.exists(file.path(data_dir, f))) "OK" else "MISSING"
      message(sprintf("  [%s] %s", status, f))
    }
    message("Optional files:")
    for (f in optional) {
      status <- if (file.exists(file.path(data_dir, f))) "OK" else "MISSING"
      message(sprintf("  [%s] %s", status, f))
    }
  }

  invisible(TRUE)
}
