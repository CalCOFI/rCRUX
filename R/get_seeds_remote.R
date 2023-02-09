#' Query primer NCBI's Blast and generate a .csv to use for blast_seeds
#'
#' @description
#' get_seeds_remote combines modified versions of [primerTree::primer_search()]
#' and primerTree's parse_primer to make [rCRUX::iterative_primer_search()]
#' which is called to query NCBI's
#' [primer BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/)
#' tool, filters the results, then aggregates them into a single data.frame.
#' It creates a directory `get_seeds_remote` in the `output_directory_path`.
#' It creates three files inside that directory. One represents the unfiltered
#' output and another represents the output after filtering with user modifiable
#' parameters and with appended taxonomy. Also generated is a summary of unique
#' taxonomic ranks after filtering.
#'
#' @details
#' get_seeds_remote passes the forward and reverse primer sequence for a given
#' PCR product to [rCRUX::iterative_primer_search()] along with the taxid(s) of
#' the organism(s) to blast, the database to search, and many additional possible
#' parameters to NCBI's primer blast tool (see Note below). Degenerate primers
#' are converted into all possible non degenerate sets and a user defined maximum
#' number of primer combinations is passed to to the API. Multiple taxids are
#' searched independently, as are multiple database searches (e.g. nt and
#' refseq_representative_genomes). The data are parsed and stored in a dataframe,
#' which is also written to a file with the suffix
#' `_unfiltered_get_seeds_remote_output.csv`.
#'
#' These hits are further filtered using [rCRUX::filter_primer_hits()] to
#' calculate and append amplicon size to the dataframe. Only hits that pass with default
#' or user modified length and number of mismatches parameters are retained.
#'
#' Taxonomy is appended to these filtered hits using
#' [rCRUX::get_taxonomizr_from_accession()]. The results are written to
#' to file with the suffix `_filtered_get_seeds_remote_output_with_taxonomy.csv`.
#' The number of unique instances for each rank in the taxonomic path for the
#' filtered hits are tallied (NAs are counted once per rank) and written to a
#' file with the suffix `_filtered_get_seeds_local_remote_taxonomic_rank_counts.txt`
#'
#'
#' Note:
#' get_seeds_remote passes many parameters to NCBI's primer blast tool.
#' You can match the parameters to the fields available in the GUI
#' [here](https://www.ncbi.nlm.nih.gov/tools/primer-blast/). First, use your
#' browser to view the page source. Search for the field you are interested in
#' by searching for the title of the field. It should be enclosed in a tag.
#' Inside the label tag, it says `for = "<name_of_parameter>"`. Copy the string
#' after for = and add it to get_seeds_remote as the name of a parameter, setting
#' it equal to whatever you like.
#'
#' As of 2022-08-16, the primer blast GUI contains some options that are not
#' implemented by [primerTree::primer_search()] and by extension [rCRUX::iterative_primer_search()]
#' primer_search doesn't include explicit documentation of allowed options, but
#' it will quickly report if an option isn't allowed, so trial and error will
#' not be very time consuming.
#'
#' Note:
#' See [rCRUX::iterative_primer_search()] and [rCRUX::modifiedPrimerTree_Functions]
#' for additional run parameters not included below.
#'
#' Check NCBI's primer blast for additional search options**
#'
#' get_seeds_remote passes many parameters to NCBI's primer blast tool. You can
#' match the parameters to the fields available in the GUI here. First, use your
#' browser to view the page source. Search for the field you are interested in
#' by searching for the title of the field. It should be enclosed in a tag.
#' Inside the label tag, it says for = "<name_of_parameter>". Copy the string
#' after for = and add it to get_seeds_remote as the name of a parameter,
#' setting it equal to whatever you like.
#'
#' As of 2022-08-16, the primer blast GUI contains some options that are not implemented by primer_search. The table below documents some of the available options.
#'
#' | Name                                   |       Default  |
#' |----------------------------------------|----------------|
#' | PRIMER_SPECIFICITY_DATABASE            | nt             |
#' | EXCLUDE_ENV                            | unchecked      |
#' | ORGANISM                               | Homo sapiens   |
#' | TOTAL_PRIMER_SPECIFICITY_MISMATCH      | 1              |
#' | PRIMER_3END_SPECIFICITY_MISMATCH       | 1              |
#' | TOTAL_MISMATCH_IGNORE                  | 6              |
#' | MAX_TARGET_SIZE                        | 4000           |
#' | HITSIZE                                | 50000          |
#' | EVALUE                                 | 30000          |
#' | WORD_SIZE                              | 7              |
#' | NUM_TARGETS_WITH_PRIMERS               | 1000           |
#' | MAX_TARGET_PER_TEMPLATE                | 100            |
#'
#'
#'
#' @param forward_primer_seq passed to [rCRUX::primer_search()], which turns it
#'        into a list of all possible non degenerate primers, then passes
#'        a user defined number of primer set combinations to NCBI.
#'        (e.g. forward_primer_seq <- "TAGAACAGGCTCCTCTAG")
#' @param reverse_primer_seq passed to [rCRUX::primer_search()], which turns it
#'        into a list of all possible non degenerate primers, then passes
#'        a user defined number of primer set combinations to NCBI.
#'        (e.g. reverse_primer_seq <-  "TTAGATACCCCACTATGC")
#' @param output_directory_path the parent directory to place the data in.
#'        (e.g. "/path/to/output/12S_V5F1_remote_111122")
#' @param metabarcode_name is passed to [rCRUX::get_seeds_remote()] which appends
#'        metabarcode_name to the beginning of each of the two files it
#'        generates (e.g. metabarcode_name <- "12S_V5F1").
#' @param accession_taxa_sql_path the path to sql created by taxonomizr
#'          (e.g. accession_taxa_sql_path <- "/my/accessionTaxa.sql")
#' @param organism a vector of character vectors. Each character vector is
#'        passed in turn to primer_search, which passes them to NCBI.
#'        get_seeds_remote aggregates all of the results into a single file.
#'        (e.g. organism = c("1476529", "7776")) - note increasing taxonomic
#'        rank (e.g. increasing from order to class) for this parameter can
#'        maximize primer hits, but can also lead to API run throttling due to
#'        memory limitations
#' @param num_permutations the number of primer permutations to search, if the
#'        degenerate bases cause more than this number of permutations to exist,
#'        this number will be sampled from all possible permutations.
#'        The default is num_permutations = 50 - Note for very degenerate bases,
#'        searches may be empty due to poor mutual matches for a given forward
#'        and reverse primer combination.
#' @param mismatch the highest acceptable mismatch value. parse_primer_hits
#'        returns a table with a mismatch column. get_seeds_remote removes each
#'        row with a mismatch greater than the specified value.
#'        The default is mismatch = 3 - Note this is smaller than [rCRUX::get_seeds_local()]
#' @param minimum_length [rCRUX::parse_primer_hits()] returns a table with a product_length
#'        column. get_seeds_remote removes each row that has a value less than
#'        minimum_length in the product_length column.
#'        The default is minimum_length = 5
#' @param maximum_length [rCRUX::parse_primer_hits()] returns a table with a
#'        product_length column. get_seeds_remote removes each row that has a
#'        value greater than maximum_length in the product_length column
#'        The default is maximum_length = 500
#' @param primer_specificity_database passed to [rCRUX::primer_search()], which passes it
#'        to NCBI.  The default is primer_specificity_database = 'nt'.
#' @param HITSIZE a primer BLAST search parameter set high to maximize the
#'        number of observations returned.
#'        The default HITSIZE = 50000 - note increasing this parameter can
#'        maximize primer hits, but can also lead to API run throttling due to
#'        memory limitations
#' @param NUM_TARGETS_WITH_PRIMERS a primer BLAST search parameter set high to
#'        maximize the number of observations returned.
#'        The default is NCBI NUM_TARGETS_WITH_PRIMERS = 1000 - - note increasing
#'        this parameter can maximize primer hits, but can also lead to API run
#'        throttling due to memory limitations
#' @param ... additional arguments passed to primer_search, see
#'        [primerTree::primer_search()] and [NCBI primer-blast tool]](https://www.ncbi.nlm.nih.gov/tools/primer-blast/)
#'        for more information.
#'
#' @return a data.frame containing the same information as the .csv it generates
#' 
#' @export
#' @examples
#'\dontrun{
#' forward_primer_seq = "TAGAACAGGCTCCTCTAG"
#' reverse_primer_seq =  "TTAGATACCCCACTATGC"
#' output_directory_path <- "/my/directory/12S_V5F1_remote_111122_modified_params"
#' metabarcode_name <- "12S_V5F1"
#' accession_taxa_sql_path <- "/my/directory/accessionTaxa.sql"
#'
#'
#' get_seeds_remote(forward_primer_seq,
#'                 reverse_primer_seq,
#'                 output_directory_path,
#'                 metabarcode_name,
#'                 accession_taxa_sql_path,
#'                 HITSIZE ='1000000',
#'                 evalue='100000',
#'                 word_size='6',
#'                 MAX_TARGET_PER_TEMPLATE = '5',
#'                 NUM_TARGETS_WITH_PRIMERS ='500000', minimum_length = 50,
#'                 MAX_TARGET_SIZE = 200,
#'                 organism = c("1476529", "7776"), return_table = FALSE)
#'
#'
#' # This results in approximately 111500 blast seed returns (there is some variation due to database updates, etc.), note the default generated approximately 1047.
#' # This assumes the user is not throttled by memory limitations.
#'}
get_seeds_remote <- function(forward_primer_seq, reverse_primer_seq,
                             output_directory_path, metabarcode_name,
                             accession_taxa_sql_path,
                             organism, mismatch = 3,
                             minimum_length = 5, maximum_length = 500,
                             primer_specificity_database = "nt", ...,
                             return_table = TRUE) {
  
  # Check paths provided
  if (!file.exists(accession_taxa_sql_path)) {
    stop("accession_taxa_sql_path does not exist.\n",
         "The path to the taxonomizr SQL file cannot be found. ",
         "Please revise the path provided:\n", accession_taxa_sql_path)
  }
  
  # Create output directories
  out <- file.path(output_directory_path, "get_seeds_local")
  dir.create(out, showWarnings = FALSE)
  
  message('Output directory: ', out, '\n')
  
  # Aggregate the primer_search return values
  # Then parse_primer_hits all of them
  raw_table <- 
    iterative_primer_search(forward_primer_seq, reverse_primer_seq,
                            organism,
                            primer_specificity_database, ...)
  
  # Throw an error if there are no results
  if (nrow(raw_table) < 1) {
    stop("Primer search returned no hits.")
  }
  
  filtered_table <- 
    filter_primer_hits(raw_table,
                       forward_primer_seq, reverse_primer_seq,
                       mismatch, minimum_length,
                       maximum_length)
  
  taxonomized_table <- 
    get_taxonomy_from_accession(filtered_table,
                                  accession_taxa_sql_path)
  
  # save output
  write.csv(taxonomized_table,
            file = file.path(out, paste0(metabarcode_name, "_filtered_get_seeds_remote_output_with_taxonomy.csv")),
            row.names = FALSE)
  
  write.csv(raw_table,
            file = file.path(out, paste0(metabarcode_name, "_unfiltered_get_seeds_remote_output.csv")),
            row.names = FALSE)
  
  # Count distinct taxonomic ranks - includes NA
  tax_rank_sum <- 
    taxonomized_table %>% 
    dplyr::summarise(
      dplyr::across(c('superkingdom', 'phylum','class','order','family','genus','species'), .fns = dplyr::n_distinct)
    )
  
  # Write output to blast_seeds_output
  write.csv(tax_rank_sum,
            file = file.path(out, paste0(metabarcode_name, "_filtered_get_seeds_remote_unique_taxonomic_rank_counts.csv")),
            row.names = FALSE)
  
  
  #return if you're supposed to
  if (return_table) {
    return(taxonomized_table)
  }
  
  invisible(NULL)
  
}
