#' Query primer_blast and generate a .csv to use for blast_seeds
#'
#' @description
#' get_seeds_remote uses a modified version of [primerTree::primer_search()] to
#' query NCBI's [primer BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/)
#' tool, filters the results, then aggregates them into a single data.frame.
#' As a side effect, it creates a directory at `output_directory_path` if one doesn't yet
#' exist, then creates a subdirectory inside `output_directory_path` named after
#' `metabarcode_name`. It creates two files inside that directory, one
#' representing the output and the other representing the output without added
#' taxonomy.
#'
#' # Additional arguments passed to primer BLAST
#'
#' get_seeds_remote passes many parameters to NCBI's primer blast tool. You can
#' match the parameters to the fields available in the GUI
#' [here](https://www.ncbi.nlm.nih.gov/tools/primer-blast/). First, use your
#' browser to view the page source. Search for the field you are interested in
#' by searching for the title of the field. It should be enclosed in a tag.
#' Inside the label tag, it says `for = "<name_of_parameter>"`. Copy the string
#' after for = and add it to get_seeds_remote as the name of a parameter, setting
#' it equal to whatever you like.
#'
#' As of 2022-08-16, the primer blast GUI
#' contains some options that are not implemented by primer_search.
#' primer_search doesn't include explicit documentation of allowed options, but
#' it will quickly report if an option isn't allowed, so trial and error will
#' not be very time consuming.
#'
#'
#' @param forward_primer_seq passed to primer_search, which turns it into a list of
#'        each primer it could be based on its degenerate primers, then passes
#'        each one in turn to NCBI
#' @param reverse_primer_seq passed to primer_search, which turns it into a list of
#'        each primer it could be based on its degenerate primers, then passes
#'        each one in turn to NCBI
#' @param output_directory_path the parent directory to place the data in.
#' @param metabarcode_name used to name the subdirectory and the files. If a
#'        directory named metabarcode_name does not exist in output_directory_path, a
#'        new directory will be created. get_seeds_remote appends
#'        metabarcode_name to the beginning of each of the two files it
#'        generates.
#' @param accession_taxa_sql_path the path to sql created by taxonomizr
#' @param organism a vector of character vectors. Each character vector is
#'        passed in turn to primer_search, which passes them to NCBI.
#'        get_seeds_remote aggregates all of the results into a single file.
#' @param mismatch the highest acceptable mismatch value. parse_primer_hits
#'        returns a table with a mismatch column. get_seeds_remote removes each
#'        row with a mismatch greater than the specified value.
#' @param minimum_length parse_primer_hits returns a table with a product_length
#'        column. get_seeds_remote removes each row that has a value less than
#'        minimum_length in the product_length column.
#' @param maximum_length parse_primer_hits returns a table with a
#'        product_length column. get_seeds_remote removes each row that has a
#'        value greater than maximum_length in the product_length column
#' @param primer_specificity_database passed to primer_search, which passes it
#'        to NCBI
#' @param HITSIZE a primer BLAST search parameter set high to maximize the
#'        number of observations returned.
#' @param NUM_TARGETS_WITH_PRIMERS a primer BLAST search parameter set high to
#'        maximize the number of observations returned.
#' @param ... additional arguments passed to primer_search, which passes it to
#'        NCBI
#' @return a data.frame containing the same information as the .csv it generates
#' @export
get_seeds_remote <- function(forward_primer_seq, reverse_primer_seq,
                            output_directory_path, metabarcode_name,
                            accession_taxa_sql_path,
                            organism, mismatch = 3,
                            minimum_length = 5, maximum_length = 500,
                            primer_specificity_database = "nt", ...,
                            return_table = TRUE) {

    # Start by making the directory and checking for the sql and whatnot.
    out <- paste0(output_directory_path, "/", metabarcode_name, "/")
    dir.create(output_directory_path)
    dir.create(out)
    if (!file.exists(accession_taxa_sql_path)) {
      stop("accession_taxa_sql_path does not exist")
    }

    # Aggregate the primer_search return values
    # Then parse_primer_hits all of them
    raw_table <- iterative_primer_search(forward_primer_seq, reverse_primer_seq,
                                          organism,
                                          primer_specificity_database, ...)
    # Throw an error if there are no results
    if (nrow(raw_table) < 1) {
      stop("Primer search returned no hits.")
    }

    filtered_table <- filter_primer_hits(raw_table,
                                          forward_primer_seq, reverse_primer_seq,
                                          mismatch, minimum_length,
                                          maximum_length)
    taxonomized_table <- get_taxonomizr_from_accession(filtered_table,
                                                        accession_taxa_sql_path)

    # save output
    save_output_as_csv(taxonomized_table,
                        "_filtered_primerTree_output_with_taxonomy", out,
                        metabarcode_name)
    save_output_as_csv(raw_table, "_raw_primerTree_output", out,
                        metabarcode_name)

    # Count distinct taxonomic ranks - includes NA
    tax_rank_sum <- dplyr::summarise_at(taxonomized_table,c('phylum','class','order','family','genus','species'),dplyr::n_distinct)

    # Write output to blast_seeds_output
    tax_rank_sum_table_path <- paste0(output_directory_path, "/", metabarcode_name, "_unique_taxonomic_rank_counts.txt")
    save_output_as_csv(tax_rank_sum, "_tax_rank_sum_table_path", out,
                        metabarcode_name)


    #return if you're supposed to
    if (return_table) {
      return(taxonomized_table)
    }
    else {
      return(NULL)
    }
}
