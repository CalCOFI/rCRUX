# We're gonna leave it up to the user whether they want to randomize the input
#' @param x a vector of fastas
#' @param sample_size the number of fastas to blast at once
#' @return a data.frame representing the output of blastn
blast_fastas <- function(x, sample_size = 1000) {
    # okay I actually need to work out the sampling on paper
    
    run_blastn()
}