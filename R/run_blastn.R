#' Runs blastn with the seed amplicon input fasta as a query
#'
#' @details
#' Calls blastn with a fasta file as the query. The user can not add
#' additional search parameters, but can modify the available parameters.
#'
#' @details
#' run_blastn takes a fasta file, and query them to a blast formatted database.
#' The result is an output table with
#' the following columns of data: accession, amplicon_length, pident,
#' query_accession, accession_sequence_length, amplicon_start, amplicon_stop,
#' sequence, evalue, BLAST_db_taxids.
#'
#' Information about the blastn parameters can be found by accessing blastn -help
#' and at [NCBI](https://www.ncbi.nlm.nih.gov/books/NBK279684/).
#'
#' @param fasta a fasta-formatted string
#' @param db a path to a directory / directories containing one or more blast-formatted database.
#'        For multiple blast databases, separate them with a space and add an extra set of quotes.
#'        (e.g blast_db_path <- "/my/ncbi_nt/nt" or blast_db_path <- '"/my/ncbi_nt/nt  /my/ncbi_ref_euk_rep_genomes/ref_euk_rep_genomes"')
#' @param temp_fasta_path a file path to write a temporary fasta to - blastn requires a file to query
#' @param ncbi_bin is the path to blast+ tools if not in the user's path.
#'        Specify only if blastn and blastdbcmd  are not in your path.
#'        The default is ncbi_bin = NULL - if not specified in path do the
#'        following: ncbi_bin = "/my/local/ncbi-blast-2.10.1+/bin/".
#' @param evalue is the number of expected hits with a similar quality score
#'        found by chance. The default is evalue = 1e-6.
#' @param coverage is the minimum percent of the query length recovered in the
#'        subject hits. The default is coverage = 50.
#' @param perID is the minimum percent identity of the query relative to the
#'        subject hits. The default is perID = 70.
#' @param align is the maximum number of subject hits to return per query
#'        blasted. The default is align = 50000.
#' @param num_threads is the number of CPUs to engage in the blastn search. The 
#'        value num_treads = 'max', uses [parallel::detectCores()] to determine 
#'        the user's number of CPUs automatically.
#' @return a tibble representing the blastn results
#' @export

run_blastn <- function(fasta, db, temp_fasta_path = NULL, ncbi_bin = NULL,
                       evalue = 1e-6, align = 50000, coverage = 50, perID = 70, num_threads = 1) {
  
  if (!is.null(ncbi_bin)){
    blastn <- file.path(ncbi_bin, 'blastn')
  } else {
    blastn = 'blastn'
  }
  
  if (num_threads == 'max') {
    cores <- parallel::detectCores()
  } else {
    cores <- num_threads
  }
  
  # This is a hacky workaround to deal with the fact
  # that blastn wants a file path as a query
  # Ideally, we would find a way (perhaps a process substitution?)
  # to pass the string directly to blastn
  # The motivation for writing the fasta here rather than elsewhere is that
  # 1) This way, the script can still generally be visualized functionally
  # 2) It allows for easy changes if we ever figure out an elegant way to do
  # the handoff
  
  if (is.null(temp_fasta_path)) {
    temp_fasta_path <- tempfile(fileext = '.fasta')
  }
  
  writeLines(fasta, con = temp_fasta_path)
  
  message("Calling blastn. This may take a long time.")
  
  blastn_output <- 
    system2(command = blastn,
            args = c("-db", db,
                     "-query", temp_fasta_path,
                     "-outfmt", paste("\"6", "saccver", "length",
                                      "pident", "qacc", "slen", "sstart",
                                      "send", "sseq", "evalue", "staxids\""),
                     "-evalue", evalue,
                     "-num_alignments", align,
                     "-qcov_hsp_perc", coverage,
                     "-perc_identity", perID,
                     "-num_threads ", cores),
            wait = TRUE,
            stdout = TRUE)
  
  
  
  file.remove(temp_fasta_path)
  
  # Format output
  column_names <- 
    c("accession",
      "amplicon_length",
      "pident",
      "query_accession",
      "accession_sequence_length",
      "amplicon_start",
      "amplicon_stop",
      "sequence",
      "evalue",
      "BLAST_db_taxids")
  
  # blastn_output is a tab-delimited string, we need split it into columns
  # as_tibble creates a one-column tibble with "value" as its col name
  blastn_output %>%
    tibble::as_tibble() %>%
    tidyr::separate(col = value, into = column_names,
                    sep = "\t")
  
}

