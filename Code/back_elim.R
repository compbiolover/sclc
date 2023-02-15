# Name: back_elim.R
# Author: Andrew Willems <awillems@vols.utk.edu>
# Purpose: Efficiently find the best predictors from a large set
# of predictors by eliminating each variable one by one based on
# their p-value
back_elim <- function(coefs = NULL,
                      return_removed_coef = TRUE,
                      run_cox = TRUE,
                      glmnet_genes = NULL,
                      data = cox_df) {
  return_list <- list()
  if (return_removed_coef == TRUE) {
    new_coefs <- filter(coefs, p_value != max(p_value) & !(genes %in% glmnet_genes))
    return_list[["Modified Coefs"]] <- new_coefs
    removed_coef <- filter(coefs, p_value == max(p_value))
    removed_coef <- removed_coef$genes
    removed_coef <- noquote(removed_coef)
    return_list[["Removed Coef"]] <- removed_coef
  }
  if (run_cox == TRUE) {
    y <- Surv(time = data$time, event = data$vital.status)
    x <- tryCatch(
      {
        as.formula(paste("y ~  ", paste(new_coefs$genes,
          collapse = "+"
        )))
      },
      error = function(e) {
        as.formula(paste("y ~  ", 1))
      }
    )
    model <- summary(coxph(data = cox_df, formula = x))
    return_list[["New Cox Model"]] <- model
  }
  return(return_list)
}
