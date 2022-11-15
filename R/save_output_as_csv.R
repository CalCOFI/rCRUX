#' Function to save data in .csv files
#'
#' @description
#' Calls write.table with file_name as the first argument
#' and `file_out_Metabarcode_description.csv` as the second
#'
#' @param file_name the object to write
#' @param description a string to include in the .csv name
#'        directly before the extension
#' @param file_out a parent directory; it should end in a /
#' @param metabarcode_name a label that appears before the file name
#' @return The value of write.table
#' @export

# This is another function that may not be necessary in the final product
save_output_as_csv <- function(file_name, description, file_out, metabarcode_name){
  write_to = paste0(file_out, metabarcode_name, description, ".csv")
  return(write.table(file_name, file = write_to, row.names=FALSE, sep = ","))
}
