#' Check BLAST+ installation
#' 
#' Checks required BLAST+ binaries are in the PATH, or found with given path. This
#' function returns a TRUE/FALSE for quick logic, and is otherwise paired with
#' [rCRUX::check_blast_plus_installation()] for detailed messaging.
#' 
#' @param ncbi_bin path to NCBI binaries if they are not in the PATH 
#' 
#' @return TRUE/FALSE with attributes ncbi_bin, checked_binaries
#' 
#' @export
has_blast_plus_binaries <- function(ncbi_bin = NULL){
  
  # A list to reference binary names
  binaries_to_check <-
    list(blastn = "blastn",
         blastdbcmd = "blastdbcmd")
  
  # Add binary paths if they are given, otherwise we expect them in the PATH
  if (!is.null(ncbi_bin)){
    binaries_to_check <- 
      lapply(binaries_to_check, function(binary){
        file.path(ncbi_bin, binary)
      })
  } 
  # -version prints 2 lines, version + build, select only version
  # return 'Not installed' if binary cannot be run
  checked_binaries <-
    c(
      blastn =
        tryCatch(system2(command = binaries_to_check$blastn, args = '-version', stdout = TRUE), 
                 error = function(e) 'Not installed')[1]
      ,
      blastdbcmd =
        tryCatch(system2(command = binaries_to_check$blastdbcmd, args = '-version', stdout = TRUE), 
                 error = function(e) 'Not installed')[1]
    )
  
  all_binaries_installed =  !all(checked_binaries %in% 'Not installed')
  
  structure(
    all_binaries_installed,
    ncbi_bin = ncbi_bin,
    checked_binaries = checked_binaries
  )

}


