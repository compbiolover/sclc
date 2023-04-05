# Name: cox_model.R
# Author: Andrew Willems <awillems@vols.utk.edu>
# Purpose: Fit an arbitrarily large number of predictors (genes) to a
# cox model to predict survival time
cox_model_fitter <- function(my_seed = 1,
                             my_alpha = 1,
                             cox_df = NULL,
                             gene_num = 1800,
                             max_it = 100000,
                             verbose = FALSE,
                             cat_preds = FALSE,
                             cox_predictors = NULL,
                             tumor_stage = FALSE,
                             tumor_n = FALSE,
                             tumor_m = FALSE,
                             n_folds = 10,
                             remove_stage = c(
                               "tumor.stage1", "tumor.stage2",
                               "tumor.stage3", "tumor.stage4"
                             ),
                             remove_n_stage = c(
                               "ajcc.n0", "ajcc.n1", "ajcc.n2",
                               "ajcc.n3"
                             ),
                             save_coefs = TRUE,
                             calc_auc = FALSE,
                             my_filename = "my_models_active_coefs.csv") {
  # Doing input sanity checks
  if (missing(my_seed)) {
    stop("You must specify a seed for reproducible analysis.")
  }

  if (missing(cox_df)) {
    stop("You must include a data.frame to carry out Cox analysis.")
  }

  if (missing(cox_predictors)) {
    stop("You must include predictors for your cox model.")
  }

  if (class(my_seed) != "numeric") {
    stop("You must specify a number for your seed.")
  }

  if (class(cox_df) != "data.frame") {
    stop("Your cox.df parameter must be of type 'data.frame'.")
  }

  # If the packages are installed, they
  # will be loaded. If they are not,
  # the packages will be installed
  # from CRAN and then loaded.
  suppressMessages(library(BiocGenerics))
  suppressMessages(library(doParallel))
  suppressMessages(library(glmnet))
  suppressMessages(library(lmtest))
  suppressMessages(library(parallel))
  suppressMessages(library(survival))
  suppressMessages(library(survminer))

  # Setting the seed for reproducible output
  set.seed(my_seed)


  # Setting the number of processors on the
  # machine to speed up the fitting
  num_of_cores <- parallel::detectCores()
  registerDoParallel(cores = num_of_cores)

  # Making the list that will store the
  # data we will return
  cox_data <- vector(mode = "list", length = 6)


  if (cat_preds == TRUE) {
    cox_predictors <- as.numeric(cox_predictors)
    my_predictors <- cox_predictors
    my_predictors <- paste("~", paste(my_predictors[1:length(my_predictors)],
      collapse = "+"
    ))
    my_x <- my_predictors
  } else {
    # The predictors for the cox model. For the SCLC we are defaulting to
    # as.character
    my_predictors <- cox_predictors
    my_predictors <- head(my_predictors, n = gene_num)
    my_predictors <- sapply(my_predictors, gsub, pattern = "-", replacement = ".")
    my_predictors <- unlist(my_predictors)
    colname_changes <- sapply(colnames(cox_df), gsub,
      pattern = "-", replacement = "."
    )
    colname_changes <- sapply(colnames(cox_df), gsub,
      pattern = "_", replacement = "."
    )
    colname_changes <- sapply(colnames(cox_df), gsub,
      pattern = "/", replacement = "."
    )
    colname_changes <- unlist(colname_changes)
    colnames(cox_df) <- colname_changes
    my_predictors <- intersect(my_predictors, colnames(cox_df))
    my_predictors <- paste("~", paste(my_predictors[1:length(my_predictors)],
      collapse = "+"
    ))
  }




  if (tumor_stage == TRUE && tumor_n == FALSE && tumor_m == FALSE) {
    my_predictors <- paste(my_predictors, "tumor.stage", sep = "+")
    my_predictors <- as.formula(my_predictors)
  } else if (tumor_stage == TRUE && tumor_n == TRUE && tumor_m == FALSE) {
    my_predictors <- paste(my_predictors, "ajcc.n", sep = "+")
    my_predictors <- as.formula(my_predictors)
  } else if (tumor_stage == TRUE && tumor_n == TRUE && tumor_m == TRUE) {
    my_predictors <- paste(my_predictors, "ajcc.m", sep = "+")
    my_predictors <- as.formula(my_predictors)
  } else {
    my_predictors <- as.formula(my_predictors)
    if (verbose == TRUE) {
      print("This is the genes only predictor....")
    }
  }

  my_x <- model.matrix(my_predictors, cox_df)


  # The response object for the cox model
  my_y <- Surv(
    time = cox_df$time,
    event = cox_df$vital.status
  )



  # The 10-fold cross-validation fit
  cv_fit <- cv.glmnet(
    x = my_x, y = my_y, nfolds = n_folds, type.measure = "C",
    maxit = max_it, family = "cox", parallel = TRUE,
    alpha = my_alpha, keep = TRUE
  )



  # Looking to see which genes are the most important
  fit <- glmnet(my_x, my_y, family = "cox", maxit = 100000)
  coefs <- coef(fit, s = cv_fit$lambda.min)



  active_index <- which(as.logical(coefs) != 0)
  active_coefficients <- coefs[active_index]
  active_genes <- rownames(coefs)[active_index]

  # Saving the coefficients of the model
  if (save_coefs == TRUE) {
    active_coefs_df <- cbind(active_genes, active_coefficients)
    write.csv(active_coefs_df, file = my_filename)
  }

  # Assessing the performance of the 10-fold cross-validation
  # on the entire data set
  if (calc_auc == TRUE) {
    model_perf <- assess.glmnet(cv_fit,
      newy = my_y,
      family = "cox", s = "lambda.min"
    )
    print(paste0(
      "The model performance on the entire data set is ",
      model_perf$C
    ))
    saveRDS(model_perf, file = "cox_model_auc.rds")
  }




  if (tumor_stage == TRUE & tumor_n == FALSE & tumor_m == FALSE) {
    print("This is the genes + tumor stage predictor...")
    active_genes <- active_genes[!active_genes %in% remove_stage]
    active_genes <- c(active_genes, "tumor.stage")
  }

  if (tumor_stage == TRUE & tumor_n == TRUE & tumor_m == FALSE) {
    print("This is the genes + tumor stage + n stage predictor...")
    active_genes <- active_genes[!active_genes %in% remove_n_stage]
    active_genes <- c(active_genes, "tumor.stage", "ajcc.n")
  }



  # Adding the relevant data bits to a list to return
  cox_data[["CV"]] <- cv_fit
  cox_data[["Coefficients"]] <- coefs
  cox_data[["Active Coefficients"]] <- active_coefficients
  cox_data[["Active Index"]] <- active_index
  cox_data[["Active Genes"]] <- active_genes
  cox_data[["Predictors"]] <- my_predictors


  # Returning our finished output
  return(cox_data)
}
