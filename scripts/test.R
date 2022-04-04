#This file will test a function for compiling the URLs from multiple primer_searches

test <- function(organism, ...) {
  url <- list()
  for(e in organism) {
    results <- primer_search(organism = e, ...)
    for(f in results) {
      url <- append(url, f$url)
    }
  }
  return(url)
}

urls <- test(c("7776", "7777"), forward = "TNGAACAGGCTCCTCTAG", reverse = "TTAGATACCCCACTATGC",
             MAX_TARGET_PER_TEMPLATE = 10)
