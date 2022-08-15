#' Function to save data in .csv files
#'
#' Calls write.table with file_name as the first argument
#' and <file_out><Metabarcode><description>.csv as the second
#'
#' @param file_name the object to write
#' @param description a string to include in the .csv name
#'        directly before the extension
#' @param file_out a parent directory; it should end in a /
#' @param Metabarcode a label that appears before the file name
#' @return The value of write.table
#' @export

# This is another function that may not be necessary in the final product
save_output_as_csv <- function(file_name, description, file_out, Metabarcode){
  write_to = paste0(file_out, Metabarcode, description, ".csv")
  return(write.table(file_name, file = write_to, row.names=FALSE, sep = ","))
}