#' Local blastn interpretation of querying primer_blast to and generate a .csv to use for blast_seeds
#'
#' @description
#' get_seeds_local is a local interpretation of [rCRUX::get_seeds_remote()] that avoids
#' querying NCBI's [primer BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/)
#' tool. Although it is slower than remotely generating blast seeds, it is not
#' subject to the arbitrary throttling of jobs that require significant memory.
#' It creates a 'get_seeds_local' directory at `output_directory_path` if one
#' doesn't yet exist, then creates a subdirectory inside `output_directory_path`
#' named after `metabarcode_name`. It creates three files inside that directory.
#' One represents the unfiltered output and another represents the output after
#' filtering with user modifiable parameters and with appended taxonomy. Also
#' generated is a summary of unique taxonomic ranks after filtering.
#'
#' @details
#' get_seeds_local passes the forward and reverse primer sequence for a given
#' PCR product to [rCRUX::run_primer_blastn()] which individually blastn queries them against a blast formatted
#' database using the task "blastn_short". The returned blast hits for each
#' sequence are matched and checked to see if they generate plausible amplicon
#' (e.g. amplify the same accession and are in the correct orientation to produce
#' a PCR product). These hits are written to a file with the suffix
#' `_unfiltered_get_seeds_local_output.csv`.  These hits are further filtered for
#' length and number of mismatches. Taxonomy is appended to these filtered hits
#' using [rCRUX::get_taxonomizr_from_accession()]. The results are written to
#' to a file with the suffix `_filtered_get_seeds_local_output_with_taxonomy.csv`.
#' The number of unique instances for each rank in the taxonomic path for the
#' filtered hits are tallied (NAs are counted once per rank) and written to a
#' file with the suffix `_filtered_get_seeds_local_unique_taxonomic_rank_counts.txt`
#'
#'
#' Note:
#' Information about the blastn parameters can be found in run_primer_blast, and
#' by accessing blastn -help. Default parameters are optimized to provide
#' results similar to that generated through remote blast via primer-blast as
#' implemented in [rCRUX::iterative_primer_search()].
#' The number of alignments returned for a given blast search is hardcoded at
#' "-num_alignments", "10000000",
#'
#' @param forward_primer_seq passed to primer_to_fasta, which turns it into fasta
#'        file to be past to get_seeds_local (e.g. forward_primer_seq <- "TAGAACAGGCTCCTCTAG").
#' @param reverse_primer_seq passed to primer_to_fasta, which turns it into fasta
#'        file to be past to get_seeds_local (e.g. reverse_primer_seq <-  "TTAGATACCCCACTATGC").
#' @param output_directory_path the parent directory to place the data in
#'        (e.g. "/path/to/output/12S_V5F1_local_111122_e300_111122").
#' @param metabarcode_name used to name the files. get_seeds_local appends
#'        metabarcode_name to the beginning of each of the files it
#'        generates (e.g. metabarcode_name <- "12S_V5F1").
#' @param accession_taxa_sql_path the path to sql created by taxonomizr
#'        (e.g. accession_taxa_sql_path <- "/my/accessionTaxa.sql")
#' @param mismatch the highest acceptable mismatch value per hit. get_seeds_local removes each
#'        row with a mismatch greater than the specified value.
#'        The default is mismatch = 6
#' @param minimum_length get_seeds_local removes each row that has a value less than
#'        minimum_length in the product_length column.
#'        The default is minimum_length = 5
#' @param maximum_length get_seeds_local removes each row that has a
#'        value greater than maximum_length in the product_length column
#'        The default is maximum_length = 500
#' @param blast_db_path path to a directory containing a blast-formatted database.
#'        (e.g blast_db_path <- "/my/ncbi_nt/nt")
#' @param task (passed to [rCRUX::run_primer_blastn()]) the task for blastn to
#'        perform - default here is "blastn_short", which is optimized
#'        for searches with queries < 50 bp
#' @param word_size (passed to [rCRUX::run_primer_blastn()]) is the fragment size
#'        used for blastn search - smaller word sizes increase sensitivity and
#'        time of the search. The default is word_size =  7
#' @param evalue (passed to [rCRUX::run_primer_blastn()]) is the number of
#'        expected hits with a similar quality score found by chance.
#'        The default is evalue = 3e-7.
#' @param coverage (passed to [rCRUX::run_primer_blastn()]) is the minimum
#'        percent of the query length recovered in the subject hits.
#'        The default is coverage = 90.
#' @param perID  (passed to [rCRUX::run_primer_blastn()]) is the minimum percent
#'        identity of the query relative to the subject hits.
#'        The default is perID = 50.
#' @param reward (passed to [rCRUX::run_primer_blastn()]) is the reward for
#'        nucleotide match. The default is reward = 2.
#' @param ncbi_bin passed to [rCRUX::run_primer_blastn()]) is the path to blast+
#'        tools if not in the user's path.  Specify only if blastn and is not in
#'        your path. The default is ncbi_bin = NULL - if not specified in path
#'        do the following: ncbi_bin = "/my/local/blast+_folder".
#' @return NULL
#' @export
#'
#' @examples
#'
#' forward_primer_seq = "TAGAACAGGCTCCTCTAG"
#' reverse_primer_seq =  "TTAGATACCCCACTATGC"
#' output_directory_path <- "/my/directory/12S_V5F1_local_111122_species_750"
#' metabarcode_name <- "12S_V5F1"
#' accession_taxa_sql_path <- "/my/directory/accessionTaxa.sql"
#' blast_db_path <- "/my/directory/ncbi_nt/nt"
#'
#'
#' get_seeds_local(forward_primer_seq,
#'                 reverse_primer_seq,
#'                 output_directory_path,
#'                 metabarcode_name,
#'                 accession_taxa_sql_path,
#'                 blast_db_path,
#'                 minimum_length = 80,
#'                 maximum_length = 150)
#'
#' # adjusting the minimum_length and maximum_length parameters reduces the number of total hits by removing reads that could result from off target amplification
#'


get_seeds_local <- function(forward_primer_seq, reverse_primer_seq,
                            output_directory_path, metabarcode_name,
                            accession_taxa_sql_path,
                            blast_db_path, mismatch = 6,
                            minimum_length = 5, maximum_length = 500,
                            primer_specificity_database = "nt", ...,
                            return_table = TRUE) {

    # Start by making the directory and checking for the sql and whatnot.
    out <- paste0(output_directory_path, "/get_seeds_local/")
    suppressWarnings(dir.create(output_directory_path))
    suppressWarnings(dir.create(out))
    if (!file.exists(accession_taxa_sql_path)) {
      stop("accession_taxa_sql_path does not exist")
    }


    # Make fasta file from primer sequences
    fasta <- paste0(">forward", "\n", forward_primer_seq, "\n", ">reverse", "\n", reverse_primer_seq, "\n")
    primer_fasta_path <- paste0(output_directory_path, "/", metabarcode_name, "primers.fasta")
    writeLines(fasta, primer_fasta_path)



    #Run run_primer on primer_fasta file
    output_table <- run_primer_blastn(primer_fasta_path, blast_db_path, ...)


    # parse amplicons from hits to forward and reverse hits
    # First isolate forward and reverse reads and rename columns
    F_only <-  dplyr::rename(dplyr::filter(output_table, qseqid=='forward'),gi=sgi, accession = saccver, mismatch_forward = mismatch, forward_start = sstart, forward_stop = send )
    R_only <-  dplyr::rename(dplyr::filter(output_table, qseqid=='reverse'), gi=sgi, accession = saccver, mismatch_reverse = mismatch, reverse_start = sstart, reverse_stop = send )

    # keep only the accessions with forward and reverse primer hits, and add a column for product length
    f_and_r <- dplyr::inner_join(x=F_only, y=R_only, by=c("accession", "gi", "staxids"))
    f_and_r <- dplyr::mutate(f_and_r, product_length=0)

    # calculate product length if F and R primer pairs are in correct orientation to make amplicon
    f_and_r <- dplyr::mutate(f_and_r, product_length = dplyr::case_when((forward_start < reverse_start & forward_start < forward_stop & reverse_stop < reverse_start ) ~ (as.numeric(reverse_start) - as.numeric(forward_start)),
                                                                    (forward_start > reverse_start & forward_start > forward_stop & reverse_stop > reverse_start) ~ (as.numeric(forward_start) - as.numeric(reverse_start)),))

    # remove all F and R primer pairs that would not make an amplicon
    f_and_r <- dplyr::filter(f_and_r, !is.na(product_length))

    #save unfiltered seeds output
    save_output_as_csv(f_and_r, "_unfiltered_get_seeds_local_output", out,
                        metabarcode_name)

    # keep only hits with acceptable product length
    f_and_r <- dplyr::filter(f_and_r, dplyr::between(product_length, minimum_length, maximum_length))

    # keep hits with accaptable number of mismatches
    f_and_r <- dplyr::filter(f_and_r, mismatch_forward <= mismatch & mismatch_reverse <= mismatch)


    taxonomized_table <- get_taxonomizr_from_accession(f_and_r,
                                                        accession_taxa_sql_path)

    # save output
    save_output_as_csv(taxonomized_table,
                        "_filtered_get_seeds_local_output_with_taxonomy", out,
                        metabarcode_name)


    # Count distinct taxonomic ranks - includes NA
    tax_rank_sum <- dplyr::summarise_at(taxonomized_table,c('superkingdom', 'phylum','class','order','family','genus','species'),dplyr::n_distinct)

    # Write output to blast_seeds_output
    tax_rank_sum_table_path <- paste0(out, "/", metabarcode_name, "_filtered_get_seeds_local_unique_taxonomic_rank_counts.txt")
    write.table(tax_rank_sum, file = tax_rank_sum_table_path, row.names = FALSE, col.names=TRUE, sep = ",")


    #return if you're supposed to
    if (return_table) {
      return(taxonomized_table)
    }
    else {
      return(NULL)
    }
}
