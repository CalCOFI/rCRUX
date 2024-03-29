#' Call primer_search with several parameters and aggregate the results
#'
#' This function acts like `primerTree::primer_search()` and
#' primerTree's parse_primer() hits all in one.
#'
#' @details
#' It calls `rCRUX::modifiedPrimerTree_Functions()` to perform tasks. Its parameters
#' are very similar to `primerTree::primer_search()`, but it takes vectors for
#' organism and for database and performs a primer search for each combination.
#' It downgrades errors from primer_search and parse_primer_hits into warnings.
#' This is useful when searching for a large number of different combinations,
#' allowing the function to output successful results.
#'
#' @param forward_primer_seq the forward primer to search (e.g. forward_primer_seq <- "TAGAACAGGCTCCTCTAG")
#' @param reverse_primer_seq the reverse primer to search (e.g. reverse_primer_seq <-  "TTAGATACCCCACTATGC")
#' @param organisms a character vector containing an id or name parseable by
#'        NCBI as an organism. If it is a vector with multiple entries, each
#'        entry will be queried separately.
#'        (e.g. organism = c("1476529", "7776"))
#' @param db which NCBI database to search. The default is db = 'nt' (e.g. nr).
#' @param ... arguments passed on to [rCRUX::primer_search()]
#' @return a data.table summarizing the results of several primer_searches
#' @export
#'

iterative_primer_search <- function(forward_primer_seq, reverse_primer_seq, organisms,
                                    db = "nt", ...) {
  output <- NULL
  # Use for loops to iterate over all the vector options
  for (org in organisms) {
    
    response <- try(
      primer_search(forward = forward_primer_seq, 
                    reverse =  reverse_primer_seq, 
                    organism = org,
                    primer_specificity_database = db, ...),
      silent = TRUE
    )
    
    
    if (inherits(response) == "try-error") {
      # To do: include useful metadata and messages
      msg <- conditionMessage(attr(response, "condition"))
      warning(msg)
    }
    else {
      # Splice the parse onto the output
      for (r in response) {
        parsed <- try(
          parse_primer_hits(r),
          silent = TRUE
        )
        if (inherits(parsed) == "try-error") {
          # To do: include useful metadata and messages
          msg <- conditionMessage(attr(response, "condition"))
          warning(msg)
          message("This occurred while processing organism ", org, ".")
        }
        else if (!is.data.frame(parsed)) {
          warning("parse_primer_hits returned an object that is not a dataframe. ", 
                  "It will be ignored.")
          message("This occurred while processing organism ", org, ".")
        }
        else {
          # turn it into a data.table
          # because I think that makes this faster?
          data.table::setDT(parsed)
          
          # It should only not be a data.frame
          # when nothing has been added yet
          # Why not initialize it as an empty data.table?
          # This is a hedge against parse_primer_hits
          # changing the output format
          if (is.null(output)) {
            output <- parsed
          }
          else {
            output <- tibble::add_row(output, parsed)
          }
        }
      }
    }
  }
  
  # Check if we got anything
  # We want it to stop because not finding anything is a problem for anything
  # that depends on it finding something and it is more helpful to simply
  # fail than return an empty data.frame that code down the line will need to
  # figure out how to deal with
  if (!is.data.frame(output)) {
    stop("Output is not a data.frame\n",
         "Hint: your searches may not have returned any results.")
  }
  
  # remove duplicate rows
  output <- dplyr::distinct(output)
  return(output)
}
