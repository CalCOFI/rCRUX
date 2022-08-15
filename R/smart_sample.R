#' Wrapper function for sample
#' If the sample size is larger than the population, returns the input vector
#' Otherwise returns the result of sampling
#' @param x a vector to sample from
#' @param size the size of the sample
#' @param ... additional arguments to pass to sample
#' @return a vector sampled from x
#' @export
smart_sample <- function(x, size, ...) {
    if (size > length(x)) {
        return(x)
    }
    else {
        return(sample(x, size, ...))
    }
}