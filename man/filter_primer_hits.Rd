% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_primer_hits.R
\name{filter_primer_hits}
\alias{filter_primer_hits}
\title{Remove accessions that do not match length and mismatch criteria determined
by parameters}
\usage{
filter_primer_hits(
  hits_table,
  forward_primer_seq,
  reverse_primer_seq,
  mismatch = 3,
  minimum_length = 5,
  maximum_length = 500
)
}
\arguments{
\item{hits_table}{a data.frame generated by parse_primer_hits
or by \code{\link[=iterative_primer_search]{iterative_primer_search()}}}

\item{forward_primer_seq}{passed to primer_search, which turns it into a list of
each primer it could be based on its degenerate primers, then passes
each one in turn to NCBI (e.g. forward_primer_seq <- "TAGAACAGGCTCCTCTAG")}

\item{reverse_primer_seq}{passed to primer_search, which turns it into a list of
each primer it could be based on its degenerate primers, then passes
each one in turn to NCBI (e.g. reverse_primer_seq <-  "TTAGATACCCCACTATGC")}

\item{mismatch}{the highest acceptable mismatch value. \code{rCRUX::parse_primer_hits()}
returns a table with a mismatch column. \code{\link[=get_seeds_remote]{get_seeds_remote()}} removes each
row with a mismatch greater than the specified value.
The default is mismatch = 3 - Note this is smaller than \code{\link[=get_seeds_local]{get_seeds_local()}}}

\item{minimum_length}{parse_primer_hits returns a table with a product_length
column. \code{\link[=get_seeds_remote]{get_seeds_remote()}} removes each row that has a value less than
minimum_length in the product_length column.
The default is minimum_length = 5}

\item{maximum_length}{parse_primer_hits returns a table with a
product_length column. \code{\link[=get_seeds_remote]{get_seeds_remote()}} removes each row that has a
value greater than maximum_length in the product_length column
The default is maximum_length = 500}
}
\value{
a data.table with problematic rows removed
}
\description{
Remove accessions that do not match length and mismatch criteria determined
by parameters
}
