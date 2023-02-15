# Name: risk_score_calculator
# Author: Andrew Willems <awillems@vols.utk.edu>
# Purpose: To calculate the risk scores for our data sets
risk_calculator <- function(my_file = active_coefs.csv,
                                  tumor_data = FALSE,
                                  n_data = FALSE,
                                  cox_df = cox_df,
                                  plot_title = "KM Plot",
                                  display_p_value = TRUE,
                                  pval_size = 5,
                                  pval_coord = c(3000, 0.9),
                                  set_test = FALSE,
                                  conf_int = TRUE,
                                  conf_int_alpha = 0.5,
                                  censor_shape = 124,
                                  censor_size = 4.5,
                                  color_pal = c("red", "blue"),
                                  risk_title = "Risk",
                                  font_title = c(40, "bold"),
                                  font_x = c(40, "bold"),
                                  font_y = c(40, "bold"),
                                  x_axis_title = "Time (days)",
                                  risk_pos = "top",
                                  risk_labs = c("High", "Low"),
                                  save_plot = TRUE,
                                  km_plot_type = "svg",
                                  plot_dpi = 600,
                                  plot_width = 7,
                                  plot_height = 7,
                                  plot_unit = "in",
                                  km_plot = "km_plot.svg") {
  # Required packages
  require(survival)

  # Custom function to modify specific aspects of 'survminer's' default theme
  custom_theme <- function() {
    theme_survminer() %+replace%
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(size = 1.5)
      )
  }

  # Since our model has coefficients that can be less than 0 we first check if
  # the median expression of the data frame is greater than 0. If it is then it
  # checks whether each gene's expression is greater than the median expression
  # of the data set it puts a value of 1 for that cell. If not the cell gets a
  # value of zero. In any other case (i.e., when the median expression of the
  # data frame is less than 0 it checks the cell and puts a value of -1 if it
  # is less than the median expression of the data frame). It puts -1 in this
  # case because the value is associated with a negative coefficient in our
  # model and we are accounting for that in our risk calculation
  risk_converter <- function(my_name = my_genes[1], my_data = risk_df,
                             my_med_exp = gene_info) {
    if (my_med_exp[counter] > 0) {
      my_data[, my_name] <- ifelse(my_data[, my_name] > my_med_exp[, my_name], 1, 0)
    } else {
      my_data[, my_name] <- ifelse(my_data[, my_name] < my_med_exp[, my_name], -1, 0)
    }
    return(my_data[, my_name])
  }

  # List to store the return objects in
  survival_return <- list()

  # Read in the data
  my_file <- read.csv(my_file)

  # Subsetting to just he columns we need
  my_file <- my_file[, 2:3]

  # Renaming them nicer names
  colnames(my_file)[c(1, 2)] <- c("gene", "coef")


  # Performing the genes only risk score calculation
  if (tumor_data == FALSE & n_data == FALSE) {
    # Simplifying the cox data frame down to just the active genes identified by
    # our model
    risk_df <- cox_df[, my_file$gene]
    risk_df <- as.matrix(risk_df)
    print(dim(risk_df))

    # Making a vector of gene signs based on coefficient values
    gene_sign <- ifelse(my_file$coef > 0, 1, -1)
    print(length(gene_sign))

    # Doing matrix multiplication of risk_df matrix by the diagonal of the
    # gene_sign vector to propagate the sign of the coefficient for each gene
    # through our risk matrix for later calculation steps
    risk_df <- risk_df %*% diag(gene_sign)
    risk_df <- as.data.frame(risk_df)

    # Giving the column names nicer names
    colnames(risk_df) <- my_file$gene

    # Making a data frame of using the apply function on the columns of the risk
    # data frame and getting their median values
    gene_info <- data.frame(med_expression = apply(risk_df, 2, median))

    # Transposing it
    gene_info <- t(gene_info)
    gene_info <- as.data.frame(gene_info)
    my_genes <- colnames(risk_df[1:length(colnames(risk_df))])

    # Getting the censoring information and survival time from the cox data frame
    risk_df$vital.status <- cox_df$vital.status
    risk_df$time <- cox_df$time

    counter <- 1
    my_converted_scores <- list()
    total_risk_length <- length(colnames(risk_df))
    risk_length <- total_risk_length - 2

    # Looping through all of the genes and comparing their expression
    for (x in my_genes) {
      current_risk <- risk_converter(
        my_name = my_genes[counter],
        my_data = risk_df[1:risk_length],
        my_med_exp = gene_info[1:risk_length]
      )
      my_converted_scores[[as.character(my_genes[counter])]] <- current_risk
      counter <- counter + 1
    }
    # Making a data frame of the converted scores by binding all the vectors
    # together
    converted_df <- data.frame(my_df = 1:dim(cox_df)[1])
    counter <- 1
    for (x in my_converted_scores) {
      converted_df <- cbind(my_converted_scores[my_genes[counter]], converted_df)
      counter <- counter + 1
    }

    # Removing the unnecessary "my_df" column from our binding of the different
    # vectors from the previous for loop
    converted_df[, "my_df"] <- NULL

    # Setting the risk score data frame's row names
    rownames(converted_df) <- rownames(cox_df)

    # Adding the survival time and censoring status
    # values from the cox data frame
    converted_df$vital.status <- cox_df$vital.status
    converted_df$time <- cox_df$time

    # Ensuring all columns in the risk data frame are of type 'numeric'
    converted_df <- apply(converted_df, 2, as.numeric)

    # Summing up each patient's risk by summing across all of the different
    # genes and their respective values for that particular patient
    patient_risks <- rowSums(x = converted_df[, 1:risk_length])
    converted_df <- as.data.frame(converted_df)
    converted_df$risk <- patient_risks

    # Ensuring all of our columns are of type 'numeric'
    converted_df <- apply(converted_df, 2, as.numeric)

    # Taking the median of the absolute value of all columns in the risk data
    # frame because we don't want the sign of the calculation to affect the
    # magnitude
    median_risk <- median(abs(converted_df[, 1:risk_length]))
    converted_df <- as.data.frame(converted_df)

    # Binarizing the risk score into high and low groups
    converted_df$risk <- ifelse(converted_df$risk > median_risk, "high", "low")

    # Fitting our inputs with the survfit function and seeing how well our risk
    # score stratifies the high versus low risk groups
    km_fit <- survfit(Surv(time, vital.status) ~ risk, data = converted_df)

    # Plotting the outcome of the fit
    finished_plot <- ggsurvplot(
      fit = km_fit,
      data = converted_df,
      pval = display_p_value,
      pval.method = set_test,
      title = plot_title,
      palette = color_pal,
      conf.int = conf_int,
      conf.int.alpha = conf_int_alpha,
      legend = risk_pos,
      legend.title = risk_title,
      legend.labs = risk_labs,
      xlab = x_axis_title,
      font.title = font_title,
      font.x = font_x,
      font.y = font_y,
      censor.shape = censor_shape,
      censor.size = censor_size,
      pval.size = pval_size,
      pval.coord = pval_coord,
      ggtheme = custom_theme()
    )

    # Viewing the finished plot
    finished_plot

    # Adding the finished plot to our return list
    survival_return[["KM Plot"]] <- finished_plot

    # Saving the plot to .svg or other specified format
    # with particular parameters
    if (save_plot == TRUE) {
      ggsave(
        filename = km_plot,
        plot = print(finished_plot$plot, newpage = FALSE),
        device = km_plot_type, dpi = plot_dpi,
        width = plot_width, height = plot_height,
        units = plot_unit
      )
    }
  }


  # Attaching the risk scoring data frame
  survival_return[["Survival DF"]] <- converted_df

  # Returning our list of finished objects
  return(survival_return)
}
