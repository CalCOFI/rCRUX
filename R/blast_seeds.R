#' Generate a metabarcode database by blasting a stratified random sample of
#' seed sequences
#'
#' @description
#' `blast_seeds()` is a wrapper function for [rCRUX::blast_datatable()] -
#' to search seeds from `get_seeds_*` functions against a blast formatted
#' database - while wrangling the output for downstream use.\cr
#'
#' It creates a permanent directory `blast_seeds_output` and
#' a temporary directory `blast_seeds_save` in the `output_directory_path`.
#' It saves from, and passes files to, [rCRUX::blast_datatable()] while the run
#' is in progress. During the final steps of the function the final data is saved in
#' `blast_seeds_output` recording the results of the blast.
#'
#' @details
#' The [rCRUX::blast_datatable()] call saves intermediate results and
#' metadata about the search as local files in the save directory generated by
#' blast_seeds. This allows the function to resume a partially
#' completed blast, mitigating the consequences of encountering an
#' error or experiencing other interruptions. To resume a partially completed
#' blast, supply the same seeds and working directory. See the documentation
#' of [rCRUX::blast_datatable()] for more information.
#'
#' During the blast_seeds the following data are cached as files and passed into
#' [rCRUX::blast_datatable()]: output_table.txt (most recent updates from the
#' blast run), blast_seeds_passed_filter.txt (seed table that tracks the blast
#' status of seeds), unsampled_indices.txt (list of seed indices that need to
#' be blasted), too_many_ns.txt (tracks seeds that have been removed due to more
#' consecutive Ns in a sequence than are acceptable (see parameter wildcards),
#' blastdbcmd_failed.txt (tracks reads that are present in the seeds database,
#' but not the local blast database. This is relevant for the results of
#' [rCRUX::get_seeds_remote()]), and lastly num_rounds.txt (tracks the number of
#' completed blast round for a given seed file).
#'
#' The final output of blast_seeds, stored in `blast_seeds_output`, are the
#' following: summary.csv (blast output with appended taxonomy),
#' {metabarcode_name}_.fasta (contains the sequence for all accessions recovered
#' during the blast search), {metabarcode_name}.taxonomy (contains the taxonomy
#' for all accessions recovered during the blast search),
#' {metabarcode_name}_blast_seeds_summary_unique_taxonomic_rank_counts.txt
#' (provides a count of all unique instances within a taxonomic rank),
#' too_many_ns.txt (tracks seeds that have been removed due to more
#' consecutive Ns in a sequence than are acceptable (see parameter wildcards),
#' blastdbcmd_failed.txt (tracks reads that are present in the seeds database,
#' but not the local blast database. This is relevant for the results of
#' [rCRUX::get_seeds_remote()]).
#'
#' @param seeds_output_path a path to a csv from `get_seeds_local()` or `get_seeds_remote()`
#'        (e.g. seeds_output_path <- '/my/rCRUX_output_directory/12S_V5F1_filtered_get_seeds_remote_output_with_taxonomy.csv')
#' @param blast_db_path a directory containing one or more blast-formatted database.
#'        For multiple blast databases, separate them with a space and add an extra set of quotes.
#'        (e.g blast_db_path <- "/my/ncbi_nt/nt" or blast_db_path <- '"/my/ncbi_nt/nt  /my/ncbi_ref_euk_rep_genomes/ref_euk_rep_genomes"')
#' @param accession_taxa_sql_path a path to the accessionTaxa sql created by
#'        taxonomizr (e.g. accession_taxa_sql_path <- "/my/accessionTaxa.sql")
#' @param output_directory_path a directory in which to save partial and complete output
#'        (e.g. "/path/to/output/12S_V5F1_local_111122_e300_111122").
#' @param metabarcode_name a prefix for the output fasta, taxonomy, and count of
#'        unique ranks.(e.g. metabarcode_name <- "12S_V5F1").
#' @param expand_vectors logical, determines whether to expand too_many_Ns
#'        and not_in db into real tables and write them in the output directory.
#'        the default is expand_vectors = TRUE.
#' @param minimum_length removes each blast result that has a value less than
#'        minimum_length in the product_length column.
#'        The default is minimum_length = 5
#' @param maximum_length removes removes each blast that has a
#'        value greater than maximum_length in the product_length column
#'        The default is maximum_length = 500
#' @param ... additional arguments passed to [rCRUX::blast_datatable()]
#'
#' @inheritDotParams blast_datatable
# blast_datatable arugments commented out
# @param sample_size passed to [rCRUX::blast_datatable()] is the the number of
#        entries to sample per rank. The default sample_size = 1 - is recommended
#        unless the user is sampling higher order taxonomy.  If there are not
#        enough seeds to sample per rank the run will end in an error.
# @param max_to_blast passed to [rCRUX::blast_datatable()] and is the maximum
#        number of entries to accumulate into a fasta before calling blastn.
#        The default is max_to_blast = 1000 - the optimal number of reads to
#        blast will depend on the user's environment (available RAM) and the
#        number of possible hits (determined by marker and parameters)
# @param wildcards passed to [rCRUX::blast_datatable()] us a character vector
#        that represents the minimum number of consecutive Ns the user will
#        tolerate in a given seed or hit sequence. The default is
#        wildcards = "NNNNNNNNNNNN"
# @param random_seed passed to [rCRUX::blast_datatable()] sets the random value generator for random stratified sampling
#        within [dplyr::slice_sample()]. Change the default (random_seed = NULL) for reproducible results.
# @param rank passed to [rCRUX::blast_datatable()] is the data column
#        representing the taxonomic rank to randomly sample. The default is
#        rank = 'genus' - sampling a lower rank  (e.g. species) will generate
#        more total hits and take more time, conversely sampling a higher rank
#        (e.g. family) will generate fewer total hits and take less time.
# @param ncbi_bin passed to [rCRUX::run_blastdbcmd()] [rCRUX::run_blastn()] is
#        the path to blast+ tools if not in the user's path.  Specify only if
#        blastn and blastdbcmd  are not in your path.
#        The default is ncbi_bin = NULL - if not specified in path do the
#        following: ncbi_bin = "/my/local/ncbi-blast-2.10.1+/bin".
# @param evalue passed to [rCRUX::run_blastn()] is the number of expected hits
#        with a similar quality score found by chance. The default is
#        evalue = 1e-6.
# @param coverage passed to [rCRUX::run_blastn()] is the minimum percent of the
#        query length recovered in the subject hits. The default is
#        coverage = 50.
# @param perID passed to [rCRUX::run_blastn()] is the minimum percent identity
#        of the query relative to the subject hits. The default is perID = 70.
# @param align passed to [rCRUX::run_blastn()] is the maximum number of subject
#        hits to return per query blasted. The default is align = 50000.
# @param minimum_length removes each row that has a value less than
#        minimum_length in the product_length column.
#        The default is minimum_length = 5
# @param maximum_length removes each row that has a
#        value greater than maximum_length in the product_length column
#        The default is maximum_length = 500
# @param num_threads number, the number of CPUs to engage in the blastn search. The
#        value 'max' can be used and which uses [parallel::detectCores()] to determine
#        the user's maximum number of CPUs automatically (use with caution; Default = 1)
#' @return NULL
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' seeds_output_path <-
#'   file.path("my/directory",
#'    "12S_V5F1_remote_111122_modified_params/blast_seeds_output",
#'    "summary.csv")
#'
#' output_directory_path <- "/my/directory/12S_V5F1_remote_111122_modified_params"
#' metabarcode_name <- "12S_V5F1"
#' accession_taxa_sql_path <- "/my/directory/accessionTaxa.sql"
#' blast_db_path <- "/my/directory/ncbi_nt/nt"
#'
#'
#' blast_seeds(seeds_output_path,
#'             blast_db_path,
#'             accession_taxa_sql_path,
#'             output_directory_path,
#'             metabarcode_name,
#'             rank = 'species',
#'             max_to_blast = 750)
#'
#' # using the rank of species will increase the number of total unique blast hits
#' # modifying the max_to_blast submits fewer reads simultaneously and reduces
#' # overall RAM while extending the run
#'}

blast_seeds <-
  function(seeds_output_path,
           blast_db_path,
           accession_taxa_sql_path,
           output_directory_path,
           metabarcode_name,
           expand_vectors = TRUE,
           minimum_length = 5,
           maximum_length = 500,
           ...) {
    dots <- list(...)
    # Setup ----
    output_dir <- file.path(output_directory_path, "blast_seeds_output")

    message('Output directory: ', output_dir, '\n')

    save_dir <- file.path(output_directory_path, "blast_seeds_save")

    # if run failed before any blast output delete save_dir
    output_table_path <- file.path(save_dir, "output_table.txt")

    if (file.exists(output_table_path) & file.size(output_table_path) < 5){
      unlink(save_dir, recursive=TRUE)
    }

    # create directories if they do not exist
    dir.create(output_directory_path, showWarnings = FALSE)
    dir.create(save_dir, showWarnings = FALSE)
    dir.create(output_dir, showWarnings = FALSE)

    # BLAST ----
    message('Blasting seeds.\n')

    blast_seeds <- utils::read.csv(seeds_output_path)

    output_table <-
      blast_datatable(blast_seeds = blast_seeds,
                      save_dir = save_dir,
                      blast_db_path = blast_db_path,
                      accession_taxa_sql_path = accession_taxa_sql_path, ...)

    message('\nBlasting complete.\n')
    message('Wrangling results.\n')

    # Wrangle ----
    # keep only hits with acceptable product length
    output_table <- dplyr::filter(output_table, dplyr::between(.data$amplicon_length, minimum_length, maximum_length))

    if (nrow(output_table) == 0){
      stop('Nothing left after filtering amplicon_length by minimum_length and maximum_length values')
    }

    # using multiple blast databases leads to duplicates so get rid of those...
    output_table <- dplyr::distinct(output_table, accession, .keep_all = TRUE)

    # Write output_table to dir/blast_seeds_output/summary.csv
    summary_csv_path <- file.path(output_dir, "summary.csv")
    utils::write.csv(output_table, file = summary_csv_path, row.names = FALSE)

    # Write a fasta
    degapped_sequence <- gsub("-", "", output_table$sequence)
    fasta <- character(length(degapped_sequence) * 2)
    fasta[c(TRUE, FALSE)] <- paste0(">", output_table$accession)
    fasta[c(FALSE, TRUE)] <- degapped_sequence

    writeLines(fasta, file.path(output_dir, paste0(metabarcode_name, ".fasta")))

    # Taxonomy file format (tidyr and dplyr)
    taxa_table <-
      output_table %>%
      dplyr::select('accession', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species') %>%
      tidyr::unite(col = 'taxonomic_path', 'superkingdom':'species', sep = ";", remove = TRUE, na.rm = FALSE) %>%
      dplyr::slice(-1)

    taxa_table_path <- file.path(output_dir, paste0(metabarcode_name, "_taxonomy.txt"))
    utils::write.table(taxa_table, file = taxa_table_path, row.names = FALSE, col.names=FALSE, sep = "\t")

    # Count distinct taxonomic ranks - includes NA
    tax_rank_sum <-
      output_table %>%
      dplyr::summarise(
        dplyr::across(.cols = c('superkingdom','phylum','class','order','family','genus','species'),
                      .fns = dplyr::n_distinct)
      )
    tax_rank_sum_table_path <- file.path(output_dir, paste0(metabarcode_name, "_blast_seeds_summary_unique_taxonomic_rank_counts.csv"))
    utils::write.csv(tax_rank_sum, file = tax_rank_sum_table_path, row.names = FALSE)


    # Read condensed vectors and expand them
    if (expand_vectors) {
      too_many_ns_path <- file.path(save_dir, "too_many_ns.txt")
      too_many_ns_indices <- as.numeric(readLines(too_many_ns_path))

      if (length(too_many_ns_indices) > 0){
        too_many_ns <- blast_seeds[too_many_ns_indices, ]
        too_many_ns_csv_path <- file.path(output_dir, "too_many_ns.csv")
        utils::write.csv(too_many_ns, file = too_many_ns_csv_path, row.names = FALSE)
      }

      blastdbcmd_failed_path <- file.path(save_dir, "blastdbcmd_failed.txt")
      blastdbcmd_failed_indices <- as.numeric(readLines(blastdbcmd_failed_path))

      if (length(blastdbcmd_failed_indices) > 0){
        blastdbcmd_failed <- blast_seeds[blastdbcmd_failed_indices, ]
        blastdbcmd_failed_csv_path <- file.path(output_dir, "blastdbcmd_failed.csv")
        utils::write.csv(blastdbcmd_failed, file = blastdbcmd_failed_csv_path, row.names = FALSE)
      }
    }

    # Clear temporary directory
    unlink(save_dir, recursive=TRUE)

    message('Done.')

    # return nothing
    invisible(NULL)
  }
