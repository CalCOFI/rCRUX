#' Check BLAST+ installation
#' 
#' Checks required BLAST+ binaries are in the PATH, or found with given path
#' 
#' @param ncbi_bin path to NCBI binaries if they are not in the PATH 
#' e.g. `/my/local/ncbi-blast-2.10.1+/bin`. Leave blank or NULL otherwise
#' @param return_success_message whether to return a message when all binaries are found
#' 
#' @return Nothing, a message or an error, depending on the outcome
check_blast_plus_installation <- function(ncbi_bin, return_success_message = FALSE){
  
  # A list to reference binary names
  binaries_to_check <-
    list(blastn = "blastn",
         blastdbcmd = "blastdbcmd")
  
  # Add binary paths if they are given, otherwise we expect them in the PATH
  if (!missing(ncbi_bin) && !is.null(ncbi_bin)){
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
  
  checked_binaries <- sub('.*: ', '', checked_binaries)
  
  checked_binaries_text <- 
    paste(names(checked_binaries), checked_binaries, sep = ': ')
  
  binaries_not_installed <- checked_binaries %in% 'Not installed'
  
  if (any(binaries_not_installed)) {
    
    if (!missing(ncbi_bin)){
      error_message <- 
        c('The NCBI binaries could not be found at ', ncbi_bin, 
          '. Please revise the path given or install NCBI Blast+')
    } else {
      error_message <-
        c('Some dependencies could not be found and require installation:\n    ',
          paste(collapse = '\n    ',
                checked_binaries_text[binaries_not_installed]),
          '\n  See README for instructions.'
        )
    }
    
    stop(error_message)
    
  }
  
  if(return_success_message){
    success_message <-
      c('Found required dependencies:\n  ',
        paste(collapse = '\n  ',
              checked_binaries_text[!binaries_not_installed])
      )
    
    message(success_message)
  }
  
  invisible(NULL)
}


