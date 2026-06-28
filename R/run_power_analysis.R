# ============================================================================
# run_power_analysis.R
#
# Runner for the a-priori feasibility analysis. Produces:
#   - Outputs/power/feasibility_report.csv  (one row per cohort/claim)
#   - Outputs/power/power_curve_survival.png
#   - Outputs/power/power_curve_correlation.png
#
# Usage (from project root):
#   source("R/run_power_analysis.R")
#   # or:  Rscript R/run_power_analysis.R
#
# This is the gate we run BEFORE any modeling: it states, given each cohort's
# size / event count, the smallest effect we could detect at 80% power. Update
# the event counts in emt_cohort_assumptions() once real clinical tables load.
# ============================================================================

source("R/power_analysis.R")

out_dir <- "Outputs/power"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Feasibility report -------------------------------------------------
report <- run_power_analysis(out_dir = out_dir)

cat("\n=========================================================\n")
cat("  A priori feasibility (alpha = 0.05, power = 0.80)\n")
cat("=========================================================\n")
print(as.data.frame(report)[, c("cohort", "claim", "n", "events",
                                "min_detectable_effect")],
      row.names = FALSE, right = FALSE)
cat(sprintf("\nReport written to %s/feasibility_report.csv\n", out_dir))

# ---- 2. Sensitivity curves (manuscript supplement) -------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  surv_grid <- power_grid_survival()
  p_surv <- plot_power_grid_survival(surv_grid)
  ggplot2::ggsave(file.path(out_dir, "power_curve_survival.png"),
                  p_surv, width = 7, height = 5, dpi = 300)

  corr_grid <- power_grid_correlation()
  p_corr <- ggplot2::ggplot(
    corr_grid,
    ggplot2::aes(x = .data$n, y = .data$power, colour = factor(.data$r))) +
    ggplot2::geom_hline(yintercept = 0.80, linetype = "dashed",
                        colour = "grey50") +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "Number of cell lines", y = "Power", colour = "|r|",
                  title = "A priori power: EMT score vs drug sensitivity",
                  subtitle = "Fisher z; two-sided alpha = 0.05") +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(file.path(out_dir, "power_curve_correlation.png"),
                  p_corr, width = 7, height = 5, dpi = 300)

  cat(sprintf("Sensitivity curves written to %s/\n", out_dir))
} else {
  cat("\n[note] Install {ggplot2} to render the sensitivity curves.\n")
}

cat("\nKey takeaways:\n")
cat("  - Single cohorts detect only HR >= ~1.5/SD; POOL cohorts for power.\n")
cat("  - Cell-line drug-response correlations are well powered for |r| >= 0.35.\n")
cat("  - The EMT x immunotherapy interaction is underpowered -> exploratory.\n")
