#' Check BLAST+ installation
#' 
#' Checks required BLAST+ binaries are in the PATH, or found with given path
#' 
#' @param ncbi_bin path to NCBI binaries if they are not in the PATH 
#' e.g. `/my/local/ncbi-blast-2.10.1+/bin`. Leave blank or NULL otherwise
#' @param return_success_message whether to return a message when all binaries are found
#' 
#' @return Nothing, a message or an error, depending on the outcome
check_blast_plus_installation <- function(ncbi_bin = NULL, return_success_message = FALSE){
  
  has_blast_plus <- has_blast_plus_binaries(ncbi_bin = ncbi_bin)
  
  ncbi_bin <- attr(has_blast_plus, 'ncbi_bin')
  checked_binaries <- sub('.*: ', '', attr(has_blast_plus, 'checked_binaries'))
  
  checked_binaries_text <- 
    paste(names(checked_binaries), checked_binaries, sep = ': ')
  
  binaries_not_installed <- checked_binaries %in% 'Not installed'

  if (any(binaries_not_installed)) {
    
    if (!is.null(ncbi_bin)){
      error_message <- 
        c('The NCBI binaries could not be found at ', ncbi_bin, 
          '. Please revise the path given or install NCBI Blast+')
    } else {
      error_message <-
        c('Some dependencies could not be found and require installation. See README for instructions.\n    ',
          paste(collapse = '\n    ',
                checked_binaries_text[binaries_not_installed])
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

