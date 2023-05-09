# Name: read_files.R
# Author: Andrew Willems <awillems@vols.utk.edu>
# Purpose: To efficiently read in all files that match a pattern in a folder and
# bind them together into a single dataframe.
read_files <- function(folder_path, pattern) {
  # Get a list of all files in the folder that match the pattern
  file_list <- list.files(path = folder_path, pattern = pattern, full.names = TRUE)
  
  # Read in all CSV files and combine into a single dataframe
  data <- do.call(rbind, lapply(file_list, read.csv))
  
  return(data)
}