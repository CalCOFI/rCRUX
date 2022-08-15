#' Query primer_blast and generate a .csv to use for rcrux_blast
#'
#' @description
#' get_blast_seeds uses a modified version of [primerTree::primer_search()] to
#' query NCBI's [primer BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/)
#' tool, filters the results, then aggregates them into a single data.frame.
#' As a side effect, it creates a directory at `file_out_dir` if one doesn't yet
#' exist, then creates a subdirectory inside `file_out_dir` named after
#' `Metabarcode_name`. It creates two files inside that directory, one
#' representing the output and the other representing the output without added
#' taxonomy.
#'
#' # Additional arguments passed to primer BLAST
#'
#' get_blast_seeds passes many parameters to NCBI's primer blast tool. You can
#' match the parameters to the fields available in the GUI here. First, use your
#' browser to view the page source. Search for the field you are interested in
#' by searching for the title of the field. It should be enclosed in a tag.
#' Inside the label tag, it says `for = "<name_of_parameter>"`. Copy the string
#' after for = and add it to get_blast_seeds as the name of a parameter, setting
#' it equal to whatever you like.
#'
#' Example: I want to set "Exon junction span" to 10. I open the source of the
#' primer designing tool and look for that string. I find the following:
#'
#' ```
#' <label class="m" for="PRIMER_ON_SPLICE_SITE">Exon junction span</label>
#' ```
#'
#' I copy PRIMER_ON_SPLICE_SITE and add it to get_blast_seeds:
#'
#' ```
#' get_blast_seeds("TAGAACAGGCTCCTCTAG", "TTAGATACCCCACTATGC",
#'              blast_seeds_parent, "12S_V5F1", accession_taxa_path,
#'              organism = c("7776"), MAX_TARGET_PER_TEMPLATE = 10,
#'              PRIMER_ON_SPLICE_SITE = "10"
#'              return_table = FALSE)
#' ```
#'
#' @param forward_primer passed to primer_search, which turns it into a list of
#'        each primer it could be based on its degenerate primers, then passes
#'        each one in turn to NCBI
#' @param reverse_primer passed to primer_search, which turns it into a list of
#'        each primer it could be based on its degenerate primers, then passes
#'        each one in turn to NCBI
#' @param file_out_dir the parent directory to place the data in.
#' @param Metabarcode_name used to name the subdirectory and the files. If a
#'        directory named Metabarcode_name does not exist in file_out_dir, a
#'        new directory will be created. get_blast_seeds appends
#'        Metabarcode_name to the beginning of each of the two files it
#'        generates.
#' @param accessionTaxa the path to sql created by taxonomizr
#' @param organism a vector of character vectors. Each character vector is
#'        passed in turn to primer_search, which passes them to NCBI.
#'        get_blast_seeds aggregates all of the results into a single file.
#' @param mismatch the highest acceptable mismatch value. parse_primer_hits
#'        returns a table with a mismatch column. get_blast_seeds removes each
#'        row with a mismatch greater than the specified value.
#' @param minimum_length parse_primer_hits returns a table with a product_length
#'        column. get_blast_seeds removes each row that has a value less than
#'        minimum_length in the product_length column.
#' @param maximum_length parse_primer_hits returns a table with a
#'        product_length column. get_blast_seeds removes each row that has a
#'        value greater than maximum_length in the product_length column
#' @param primer_specificity_database passed to primer_search, which passes it
#'        to NCBI
#' @param ... additional arguments passed to primer_search, which passes it to
#'        NCBI
#' @return a data.frame containing the same information as the .csv it generates
#' @export
get_blast_seeds <- function(forward_primer, reverse_primer,
                            file_out_dir, Metabarcode_name,
                            accessionTaxa,
                            organism, mismatch = 3,
                            minimum_length = 5, maximum_length = 500,
                            primer_specificity_database = "nt", ...,
                            return_table = TRUE) {

    # Start by making the directory and checking for the sql and whatnot.
    out <- paste0(file_out_dir, "/", Metabarcode_name, "/")
    dir.create(file_out_dir)
    dir.create(out)
    if (!file.exists(accessionTaxa)) {
      stop("accessionTaxa does not exist")
    }

    # Aggregate the primer_search return values
    # Then parse_primer_hits all of them
    raw_table <- iterative_primer_search(forward_primer, reverse_primer,
                                          organism,
                                          primer_specificity_database, ...)
    # Throw an error if there are no results
    if (nrow(raw_table) < 1) {
      stop("Primer search returned no hits.")
    }

    filtered_table <- filter_primer_hits(raw_table,
                                          forward_primer, reverse_primer,
                                          mismatch, minimum_length,
                                          maximum_length)
    taxonomized_table <- get_taxonomizr_from_accession(filtered_table,
                                                        accessionTaxa)

    # save output
    save_output_as_csv(taxonomized_table,
                        "_primerTree_output_with_taxonomy", out,
                        Metabarcode_name)
    save_output_as_csv(raw_table, "_raw_primerTree_output", out,
                        Metabarcode_name)

    #return if you're supposed to
    if (return_table) {
      return(taxonomized_table)
    }
    else {
      return(NULL)
    }
}