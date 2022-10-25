#' Runs blastn with the input primer sequences converted to a primer_fasta file
#' as a query
#'
#' Takes the input string, writes it to a temporary file,
#' and calls blastn with that file as the query. Warning:
#' if a file at the path specified by temp already exists,
#' it will be overwritten then deleted.
#'
#' @param primer_fasta path to the primer fasta file
#' @param ncbi_bin if not null use it as the parent directory for blastn
#' @param db path to blast formatted database
#' @param task the task for blastn to perform - default here is "blastn_short",
#' which is optimized for searches with queries < 50 bp
#' @param word_size is the fragment size used for blastn search - smaller word
#' sizes increase sensitivity and time of the search - default value is 7
#' @param evalue is the number of expected hits with a similar quality score
#' found by chance - default is 3e-7.
#' @param coverage is the minimum percent of the query length recovered in the
#' subject hits
#' @param perID is the minimum percent identity of the query relative to the
#' subject hits, the default is 50
#' @param reward is the reward for nucleotide match, the default is 2
#' @return a tibble representing the blastn results
#' @export

#might have to hardcode align...  see if the package works diff than terminal

run_primer_blastn <- function(primer_fasta, db, ncbi_bin = NULL, task = "blastn-short", word_size = 7,
                       evalue = '3e+07', coverage = 90, perID = 50, reward = 2) {


  # Determine arguments
  cores <- parallel::detectCores()

  message("Calling blastn for primers. This may take a long time.")

  # System call

  if (is.null(ncbi_bin)) {
    blastn_output <- system2(command = "blastn",
                             args = c("-db", db,
                                      "-task", task,
                                      "-query", primer_fasta,
                                      "-outfmt", paste("\"6", "qseqid", "sgi",
                                                       "saccver", "mismatch", "sstart",
                                                       "send", "staxids\""),
                                      "-evalue", evalue,
                                      "-num_alignments", 10000000,
                                      "-qcov_hsp_perc", coverage,
                                      "-perc_identity", perID,
                                      "-reward", reward,
                                      "-word_size", word_size,
                                      "-num_threads ", cores),
                             wait = TRUE,
                             stdout = TRUE)
  } else {
    blastn <- paste0(ncbi_bin, "blastn")
    blastn_output <- system2(command = "blastn",
                             args = c("-db", db,
                                      "-task", task,
                                      "-query", primer_fasta,
                                      "-outfmt", paste("\"6", "qseqid", "sgi",
                                                       "saccver", "mismatch", "sstart",
                                                       "send", "staxids\""),
                                      "-evalue", evalue,
                                      "-num_alignments", align,
                                      "-qcov_hsp_perc", coverage,
                                      "-perc_identity", perID,
                                      "-reward", reward,
                                      "-word_size", word_size,
                                      "-num_threads ", cores),
                             wait = TRUE,
                             stdout = TRUE)
  }

  file.remove(primer_fasta)

  # Format output
  column_names <-  c("qseqid",
                     "sgi",
                     "saccver",
                     "mismatch",
                     "sstart",
                     "send",
                     "staxids")


  output_table <- blastn_output %>%
    tibble::as_tibble() %>%
    tidyr::separate(col = value, into = column_names,
                    sep = "\t")
  return(output_table)

}




`%>%` <- magrittr::`%>%`
