#' Run blastn with a fasta file
#'
#' @details
#' Calls blastn with a primer fasta file as the query. The user can not add
#' additional search parameters, but can modify the available parameters.
#'
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
#'
#' @param primer_fasta path to the primer fasta file
#' @param db a path to a directory / or directories containing one or more blast-formatted database.
#'        For multiple blast databases, separate them with a space and add an extra set of quotes.
#'        (e.g blast_db_path <- "/my/ncbi_nt/nt" or blast_db_path <- '"/my/ncbi_nt/nt  /my/ncbi_ref_euk_rep_genomes/ref_euk_rep_genomes"')
#' @param task the task for blastn to perform. The default here is "blastn-short",
#'        which is optimized for searches with queries < 50 bp.
#' @param word_size is the fragment size used for blastn search. Smaller word
#'        sizes increase sensitivity and time of the search. The default is word_size =  7.
#' @param evalue is the number of expected hits with a similar quality score
#'        found by chance. The default is evalue = '3e-7'.
#' @param coverage is the minimum percent of the query length recovered in the
#'        subject hits. The default is coverage = 90.
#' @param perID is the minimum percent identity of the query relative to the
#'        subject hits. The default is perID = 50.
#' @param align is the maximum number of subject hits to return per query
#'        blasted. The default is align = '10000000'. - to few alignments will
#'        result in no matching pairs of forward and reverse primers.  To many
#'        alignments can result in an error due to RAM limitations.
#' @param num_threads number, the number of CPUs to engage in the blastn search. The
#'        value 'max' can be used and which uses [parallel::detectCores()] to determine
#'        the user's maximum number of CPUs automatically (use with caution; Default = 1)
#' @param reward is the reward for nucleotide match. The default is reward = 2.
#' @param ncbi_bin is the path to blast+ tools if not in the user's PATH
#'        Specify only if blastn and blastdbcmd  are not in your path.
#'        The default is ncbi_bin = NULL - if not specified in path do the
#'        following: ncbi_bin = "/my/local/ncbi-blast-2.10.1+/bin/".
#'
#'
#' @return a tibble 'output_table' representing the blastn results
#'
#' @examples
#'
#' temp_fasta <- tempfile(fileext = '.fasta')
#'
#' test_primers <-
#'   c('>primer_forward', 'AGAGGAGCGCGGAATTCC',
#'   '>primer_reverse', 'TACCTTGTTACGACTT')
#'
#' writeLines(test_primers, temp_fasta)
#'
#' blast_db_path <- file.path(system.file(package = 'rCRUX', 'mock-db/blastdb'), 'mock-db')
#'
#' # Returns a data.frame of results
#' result <- run_primer_blastn(primer_fasta = temp_fasta, db = blast_db_path)
#'
#' @export
run_primer_blastn <-
  function(primer_fasta,
           db,
           task = "blastn-short",
           word_size = 7,
           evalue = '3e+07',
           align = '10000000',
           coverage = 90,
           perID = 50,
           reward = 2,
           num_threads = 1,
           ncbi_bin = NULL) {

    # Assume paths to things have been checked?
    check_blast_plus_installation(ncbi_bin = ncbi_bin)
    check_blast_db(db)

    # Not sure this is the best default
    if (num_threads == 'max') {
      cores <- parallel::detectCores()
    } else {
      cores <- num_threads
    }

    # Prepare call to blastn
    if (!is.null(ncbi_bin)){
      blastn <- file.path(ncbi_bin, 'blastn')
    } else {
      blastn = 'blastn'
    }

    args <-
      c("-db", db,
        "-task", task,
        "-query", primer_fasta,
        "-outfmt", '"6 qseqid sgi saccver mismatch sstart send staxids"',
        "-evalue", evalue,
        "-num_alignments", align,
        "-qcov_hsp_perc", coverage,
        "-perc_identity", perID,
        "-reward", reward,
        "-word_size", word_size,
        "-num_threads", cores,
        "-mt_mode", 1)

    message("Calling blastn for primers. This may take a long time.\n")
    message(paste(blastn, paste(args, collapse = ' ')))

    # Catch stdout to character vector (a tab-delimited table)
    blastn_output <-
      system2(command = blastn,
              args = args,
              wait = TRUE,
              stdout = TRUE)

    file.remove(primer_fasta)

    # Wrangle return data to a tibble
    column_names <-
      c("qseqid",
        "sgi",
        "saccver",
        "mismatch",
        "sstart",
        "send",
        "staxids")


    blastn_output %>%
      tibble::as_tibble() %>%
      tidyr::separate(col = .data$value,
                      into = column_names,
                      sep = "\t")

  }
