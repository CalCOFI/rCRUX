#' Runs blastn with the input primer sequences converted to a primer_fasta file
#' as a query
#'
#' @details
#' Calls blastn with a primer fasta file as the query. The user can not add
#' additional search parameters, but can modify the available parameters.
#'
#' @details
#' run_primer_blastn takes a fasta file containing primers, uses blastn-short to
#' query them to a blast formatted database. The result is an output table with
#' the following columns of data: qseqid (query subject id), sgi (subject gi),
#' saccver (subject accession version), mismatch (number of mismatches between
#' the subject a query), sstart (subject start), send (subject end), staxids
#' (subject taxids).
#'
#' Information about the blastn parameters can be found by accessing blastn -help
#' and at [NCBI](https://www.ncbi.nlm.nih.gov/books/NBK279684/).
#'
#' Note:
#' The number of alignments returned for a given blast search is hardcoded at
#' "-num_alignments 10000000".
#'
#' @param primer_fasta path to the primer fasta file
#' @param ncbi_bin is the path to blast+ tools if not in the user's path.
#'        Specify only if blastn and blastdbcmd  are not in your path.
#'        The default is ncbi_bin = NULL - if not specified in path do the
#'        following: ncbi_bin = "/my/local/blast+_folder".
#' @param db path to blast formatted database
#' @param task the task for blastn to perform - default here is "blastn_short",
#'        which is optimized for searches with queries < 50 bp
#' @param word_size is the fragment size used for blastn search - smaller word
#'        sizes increase sensitivity and time of the search.
#'        The default is word_size =  7
#' @param evalue is the number of expected hits with a similar quality score
#'        found by chance. The default is evalue = 3e-7.
#' @param coverage is the minimum percent of the query length recovered in the
#'        subject hits. The default is coverage = 90.
#' @param perID is the minimum percent identity of the query relative to the
#'        subject hits. The default is perID = 50.
#' @param reward is the reward for nucleotide match. The default is reward = 2.
#' @return a tibble 'output_table' representing the blastn results
#' @export


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
                                      "-num_alignments", "10000000",
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
                                      "-num_alignments", "10000000",
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
