#' Local blastn interpretation of querying primer_blast to and generate a .csv to use for rcrux_blast
#'
#' @description
#' rCRUX_primer_blast is a local interpretation of get_blast_seeds that avoids
#' querying NCBI's [primer BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/)
#' tool. Although it is slower than remotly genertating blast seeds, it is not
#' subject to the arbitrary throttling of jobs that require significant memory.
#' It creates a directory at `file_out_dir` if one doesn't yet
#' exist, then creates a subdirectory inside `file_out_dir` named after
#' `Metabarcode_name`. It creates two files inside that directory, one
#' representing the output and the other representing the output without added
#' taxonomy.
#'
#' Information about the blastn parameters can be found in run_primer_blast, and
#' by accessing blastn -help.  Default parameters were optimized to provide
#' results similar or expanded results that through remote blast via primer-blast.
#'
#' @param forward_primer passed to primer_to_fasta, which turns it into fasta
#'        file to be past to rCRUX_primer_blast
#' @param reverse_primer passed to primer_to_fasta, which turns it into fasta
#'        file to be past to rCRUX_primer_blast
#' @param file_out_dir the parent directory to place the data in.
#' @param Metabarcode_name used to name the subdirectory and the files. If a
#'        directory named Metabarcode_name does not exist in file_out_dir, a
#'        new directory will be created. get_blast_seeds appends
#'        Metabarcode_name to the beginning of each of the two files it
#'        generates.
#' @param accessionTaxa the path to sql created by taxonomizr
#' @param mismatch the highest acceptable mismatch value per hit. rCRUX_primer_blast removes each
#'        row with a mismatch greater than the specified value.
#' @param minimum_length rCRUX_primer_blast removes each row that has a value less than
#'        minimum_length in the product_length column.
#' @param maximum_length rCRUX_primer_blast removes each row that has a
#'        value greater than maximum_length in the product_length column
#' @param db a directory with a blast-formatted database
#' @param task the task for blastn to perform - default here is "blastn_short",
#'        which is optimized for searches with queries < 50 bp
#' @param word_size is the fragment size used for blastn search - smaller word
#'        sizes increase sensitivity and time of the search - default value is 7
#' @param evalue is the number of expected hits with a similar quality score
#'        found by chance - default is 3e-7.
#' @param coverage is the minimum percent of the query length recovered in the
#'        subject hits
#' @param perID is the minimum percent identity of the query relative to the
#'        subject hits, the default is 50
#' @param reward is the reward for nucleotide match, the default is 2
#' @return a data.frame containing the same information as the .csv it generates
#' @export


rcrux_primer_blast <- function(forward_primer, reverse_primer,
                            file_out_dir, Metabarcode_name,
                            accessionTaxa,
                            db, mismatch = 6,
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


    # Make fasta file from primer sequences
    fasta <- paste0(">forward", "\n", forward, "\n", ">reverse", "\n", reverse, "\n")
    primer_fasta_path <- paste0(file_out_dir, "/", Metabarcode_name, "primers.fasta")
    writeLines(fasta, primer_fasta_path)



    #Run run_primer on primer_fasta file
    output_table <- run_primer_blastn(primer_fasta_path, db, ...)


    # parse amplicons from hits to forward and reverse hits
    # First isolate forward and reverse reads and rename columns
    F_only <-  dplyr::rename(dplyr::filter(output_table, qseqid=='forward'),gi=sgi, accession = saccver, mismatch_forward = mismatch, forward_start = sstart, forward_end = send )
    R_only <-  dplyr::rename(dplyr::filter(output_table, qseqid=='reverse'), gi=sgi, accession = saccver, mismatch_reverse = mismatch, reverse_start = sstart, reverse_end = send )

    # keep only the accessions with forward and reverse primer hits, and add a column for product length
    f_and_r <- dplyr::inner_join(x=F_only, y=R_only, by=c("accession", "gi", "staxids"))
    f_and_r <- dplyr::mutate(f_and_r, product_length=0)

    # calculate product length if F and R primer pairs are in correct orientation to make amplicon
    f_and_r <- dplyr::mutate(f_and_r, product_length = dplyr::case_when((forward_start < reverse_start & forward_start < forward_end & reverse_end < reverse_start ) ~ (as.numeric(reverse_start) - as.numeric(forward_start)),
                                                                    (forward_start > reverse_start & forward_start > forward_end & reverse_end > reverse_start) ~ (as.numeric(forward_start) - as.numeric(reverse_start)),))

    # remove all F and R primer pairs that would not make an amplicon
    f_and_r <- dplyr::filter(f_and_r, !is.na(product_length))

    # keep only hits with acceptable product length
    f_and_r <- dplyr::filter(f_and_r, dplyr::between(product_length, minimum_length, maximum_length))

    # keep hits with accaptable number of mismatches
    f_and_r <- dplyr::filter(f_and_r, mismatch_forward <= mismatch & mismatch_reverse <= mismatch)


    taxonomized_table <- get_taxonomizr_from_accession(f_and_r,
                                                        accessionTaxa)

    # save output
    save_output_as_csv(taxonomized_table,
                        "_primerTree_output_with_taxonomy", out,
                        Metabarcode_name)
    save_output_as_csv(f_and_r, "_rcrux_primer_blast_output", out,
                        Metabarcode_name)

    # Count distinct taxonomic ranks - includes NA
    tax_rank_sum <- dplyr::summarise_at(taxonomized_table,c('phylum','class','order','family','genus','species'),dplyr::n_distinct)

    # Write output to rcrux_blast_output
    tax_rank_sum_table_path <- paste0(file_out_dir, "/", Metabarcode_name, "_unique_taxonomic_rank_counts.txt")
    save_output_as_csv(tax_rank_sum, "_tax_rank_sum_table_path", out,
                        Metabarcode_name)


    #return if you're supposed to
    if (return_table) {
      return(taxonomized_table)
    }
    else {
      return(NULL)
    }
}
