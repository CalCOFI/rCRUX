#' Local blastn interpretation of querying primer_blast to and generate a .csv to use for blast_seeds
#'
#' @description
#' get_seeds_local is a local interpretation of [rCRUX::get_seeds_remote()] that avoids
#' querying NCBI's [primer BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/)
#' tool. Although it is slower than remotely generating blast seeds, it is not
#' subject to the arbitrary throttling of jobs that require significant memory.
#' It creates a 'get_seeds_local' directory at `output_directory_path` if one
#' doesn't yet exist, then creates a subdirectory inside `output_directory_path`
#' named after `metabarcode_name`. It creates four permenant files inside that
#' directory. Several temporary files are generated during the run for keeping
#' track of primers being blasted and needing to be blasted, and raw blast output.
#' One represents the unfiltered output and another represents the output after
#' filtering with user modifiable parameters and with appended taxonomy. Also
#' generated is a summary of unique taxonomic ranks after filtering, and a
#' fasta of the primers used for blast.
#'
#' @details
#' get_seeds_local passes the forward and reverse primer sequence for a given
#' PCR product to [rCRUX::run_primer_blastn()]. In the case of a non degenerate
#' primer set only two primers will be passed to run_primer_blast.  In the case
#' of a degenerate primer set, get_seeds_local will get all possible versions of
#' the degenerate primer(s) (using primerTree's enumerate_primers() function),
#' randomly sample a user defined number of forward and reverse primers, and
#' generate a fasta file. The selected primers are subset and passed to
#' run_primer_blastn which queries each primer against a blast formatted database
#' using the task "blastn_short".
#'
#' Output is cashed after each sucessful run of run_primer_blastn, so if a run
#' is interrupted the user can resubmit the command and pick up where they left
#' off.  The user can modify parameters for the run with the exception of
#' num_fprimers_to_blast and num_rprimers_to_blast.
#'
#' The returned blast hits for each
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
#'
#' @param forward_primer_seq which which turns degenerate primers into into a
#'        list of all possible non degenerate primers and converts the primer(s)
#'        into to a fasta file to be past to run_primer_blastn.
#'        (e.g. forward_primer_seq <- "TAGAACAGGCTCCTCTAG" or
#'        forward_primer_seq <- c("TAGAACAGGCTCCTCTAG", "GGWACWGGWTGAACWGTWTAYCCYCC")
#' @param reverse_primer_seq which which turns degenerate primers into into a
#'        list of all possible non degenerate primers and converts the primer(s)
#'        into to a fasta file to be past to run_primer_blastn.
#'        (e.g reverse_primer_seq <-  "TTAGATACCCCACTATGC" or
#'        reverse_primer_seq <- c("TTAGATACCCCACTATGC", "TANACYTCNGGRTGNCCRAARAAYCA")
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
#' @param blast_db_path a directory containing one or more blast-formatted database.
#'        For multiple blast databases, separate them with a space and add an extra set of quotes.
#'        (e.g blast_db_path <- "/my/ncbi_nt/nt" or blast_db_path <- '"/my/ncbi_nt/nt  /my/ncbi_ref_euk_rep_genomes/ref_euk_rep_genomes"')
#' @param task (passed to [rCRUX::run_primer_blastn()]) the task for blastn to
#'        perform - default here is "blastn_short", which is optimized
#'        for searches with queries < 50 bp
#' @param word_size (passed to [rCRUX::run_primer_blastn()]) is the fragment size
#'        used for blastn search - smaller word sizes increase sensitivity and
#'        time of the search. The default is word_size =  7
#' @param evalue (passed to [rCRUX::run_primer_blastn()]) is the number of
#'        expected hits with a similar quality score found by chance.
#'        The default is evalue = '3e-7'.
#' @param coverage (passed to [rCRUX::run_primer_blastn()]) is the minimum
#'        percent of the query length recovered in the subject hits.
#'        The default is coverage = 90.
#' @param perID  (passed to [rCRUX::run_primer_blastn()]) is the minimum percent
#'        identity of the query relative to the subject hits.
#'        The default is perID = 50.
#' @param reward (passed to [rCRUX::run_primer_blastn()]) is the reward for
#'        nucleotide match. The default is reward = 2.
#' @param align is the maximum number of subject hits to return per query
#'        blasted. The default is align = '10000000'. - to few alignments will
#'        result in no matching pairs of forward and reverse primers.  To many
#'        alignments can result in an error due to RAM limitations.
#' @param num_fprimers_to_blast is the maximum number of possible forward primers
#'        to blast. This is relevant for degenerate primers, all possible primers
#'        from a degenerate sequence are enumerated, and the user can choose a
#'        number to be randomly sampled and used for primer blast.
#'        The default is num_fprimers_to_blast = 50
#' @param num_rprimers_to_blast is the maximum number of possible reverse primers
#'        to blast. This is relevant for degenerate primers, all possible primers
#'        from a degenerate sequence are enumerated, and the user can choose a
#'        number to be randomly sampled and used for primer blast.
#'        The default is num_rprimers_to_blast = 50
#' @param max_to_blast is the number of primers to blast simultaneously.
#'        The default is max_to_blast = 2. - Increasing this number will decrease
#'        overall run time, but increase the amount of RAM required.
#' @param ncbi_bin passed to [rCRUX::run_primer_blastn()]) is the path to blast+
#'        tools if not in the user's path.  Specify only if blastn and is not in
#'        your path. The default is ncbi_bin = NULL - if not specified in path
#'        do the following: ncbi_bin = "/my/local/ncbi-blast-2.10.1+/bin".
#' @return NULL
#' @export
#'
#' @examples
#'
#' # Non degenerate primer example: 12S_V5F1 (Riaz et al. 2011)
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
#'
#' # Degenerate primer example - mlCOIintF/jgHC02198 (Leray et al. 2013)
#' # Note: this will take considerable time and computational resources
#'
#' forward_primer_seq <- "GGWACWGGWTGAACWGTWTAYCCYCC"
#' reverse_primer_seq <- "TANACYTCNGGRTGNCCRAARAAYCA"
#' output_directory_path <- "/my/directory/CO1_local"
#' metabarcode_name <- "CO1"
#'
#'
#' get_seeds_local(forward_primer_seq,
#'                 reverse_primer_seq,
#'                 output_directory_path,
#'                 metabarcode_name,
#'                 accession_taxa_sql_path,
#'                 blast_db_path,
#'                 minimum_length = 200,
#'                 maximum_length = 400,
#'                 aligns = '10000'
#'                 num_rprimers_to_blast = 200,
#'                 num_rprimers_to_blast = 2000,
#'                 max_to_blast = 10)
#'
#'
#' # Non Degenerate but high return primer example - 18S (Amaral-Zettler et al. 2009)
#' # Note: this will take considerable time and computational resources
#'
#' forward_primer_seq <- "GTACACACCGCCCGTC"
#' reverse_primer_seq <- "TGATCCTTCTGCAGGTTCACCTAC"
#' output_directory_path <- "/my/directory/18S_local"
#' metabarcode_name <- "18S"
#'
#'
#' get_seeds_local(forward_primer_seq,
#'                 reverse_primer_seq,
#'                 output_directory_path,
#'                 metabarcode_name,
#'                 accession_taxa_sql_path,
#'                 blast_db_path,
#'                 minimum_length = 250,
#'                 maximum_length = 350,
#'                 max_to_blast = 1)
#'
#' # blasting two primers at a time can max out a system's RAM, however blasting one at a time is more feasable for personal computers with 16 GB RAM
#'
#'
#'





get_seeds_local <- function(forward_primer_seq, reverse_primer_seq,
                            output_directory_path, metabarcode_name,
                            accession_taxa_sql_path,
                            blast_db_path, mismatch = 6,
                            minimum_length = 5, maximum_length = 500,
                            primer_specificity_database = "nt",
                            num_fprimers_to_blast = 50,
                            num_rprimers_to_blast = 50, max_to_blast = 2,
                            align = '10000000',  ...,
                            return_table = TRUE) {

  # Start by making the directory and checking for the sql and whatnot.
  out <- paste0(output_directory_path, "/get_seeds_local/")
  suppressWarnings(dir.create(output_directory_path))
  suppressWarnings(dir.create(out))
  if (!file.exists(accession_taxa_sql_path)) {
    stop("accession_taxa_sql_path does not exist")
  }

  # make fasta file paths for storing primer blast runs
  fasta_path <- paste0(out, metabarcode_name, "_primers_selected_for_blastn.fasta")
  to_blast_path <- paste0(out, metabarcode_name, "_subset_for_blastn.fasta")
  left_to_blast_path <- paste0(out, metabarcode_name, "_need_to_blastn.fasta")
  append_table_path <- paste0(out, metabarcode_name, "_temp_blast_output.csv")

  # with any luck this function will pick up where it left off by checking for the left_to_blast file.
  # If it does not exist it starts from scratch.  If it does not exist it skips the beginning bit and go strait to sub-setting for blast.

  if (!file.exists(left_to_blast_path)) {

    #key for degenerate bases, copied from modifiedPrimerTreeFunctions
    iupac = list( "M" = list("A", "C"),
                  "R" = list("A", "G"),
                  "W" = list("A", "T"),
                  "S" = list("C", "G"),
                  "Y" = list("C", "T"),
                  "K" = list("G", "T"),
                  "V" = list("A", "C", "G"),
                  "H" = list("A", "C", "T"),
                  "D" = list("A", "G", "T"),
                  "B" = list("C", "G", "T"),
                  "N" = list("A", "C", "G", "T"),
                  "I" = list("A", "T", "C"))



    # add multiple primers and or get the possible combinations for degenerate primers using functions in primer tree
    fPrimer <- NULL
    rPrimer <- NULL

    for (pf in forward_primer_seq ) {
      forward_primers = enumerate_ambiguity(pf)
      forward.df <- data.frame(forward=forward_primers, stringsAsFactors = FALSE)
      fPrimer <- rbind(fPrimer, forward.df)
    }

    for (pr in reverse_primer_seq ) {
      reverse_primers = enumerate_ambiguity(pr)
      reverse.df <- data.frame(reverse=reverse_primers, stringsAsFactors = FALSE)
      rPrimer <- rbind(rPrimer, reverse.df)
    }


    # count rows / number of primers
    nforward <- nrow(fPrimer)
    nreverse <- nrow(rPrimer)


    #subset primers if user so chooses
    if(nforward > num_fprimers_to_blast){
      message("")
      message(paste0( 'Forward primers have ', nforward, ' possible sequences due to degenerate bases. Randomly sampling ', num_fprimers_to_blast, ' forward primers. To change this, modify num_fprimers_to_blast.'  ))
      forward_sample = dplyr::sample_n(fPrimer, num_fprimers_to_blast, replace=FALSE)
    } else {
      message("")
      message(paste0(  nforward, ' forward primer(s) will be blasted.'))
      forward_sample = fPrimer
    }

    if(nreverse > num_rprimers_to_blast){
      message(paste0( 'Reverse primers have ', nreverse, ' possible sequences due to degenerate bases. Randomly sampling ', num_rprimers_to_blast, ' reverse primers. To change this, modify num_rprimers_to_blast.'  ))
      message("")
      reverse_sample = dplyr::sample_n(rPrimer, num_rprimers_to_blast, replace=FALSE)
    }  else {
      message(paste0( nreverse, ' reverse primer(s) will be blasted.'  ))
      message("")
      reverse_sample = rPrimer
    }


    # forward fasta first
    fastaf <- character(nrow(forward_sample) * 2 )
    fastaf[c(TRUE, FALSE)] <- paste0(">forward_row_", row.names(forward_sample))
    fastaf[c(FALSE, TRUE)] <- forward_sample$forward

    # then reverse fasta
    fastar <- character(nrow(reverse_sample) *2 )
    fastar[c(TRUE, FALSE)] <- paste0(">reverse_row_", row.names(reverse_sample))
    fastar[c(FALSE, TRUE)] <- reverse_sample$reverse

    # write forward to fasta output then append reverse
    writeLines(fastaf, fasta_path)
    write(fastar, file = fasta_path, append = TRUE)

    # column names for the blast output
    column_names <-  c("qseqid",
                       "sgi",
                       "saccver",
                       "mismatch",
                       "sstart",
                       "send",
                       "staxids")

    # make a tibble to store blast output and add column names
    append_table <- tibble::as_tibble(data.frame(matrix(nrow=0,ncol=length(column_names))))
    colnames(append_table) <- column_names


  } else {

    fasta_path = left_to_blast_path
    append_table = read.csv(append_table_path, colClasses = "character")

  }

  #
  input <- readr::read_lines(fasta_path)


  message(paste0('Reads will be blasted in subsets of up to ', max_to_blast,  ' read(s). To change this, modify max_to_blast.'))
  message("")

  # take a subset of the primers and blast - keep subsetting until we finish
  while (length(input) > 0){
    # lines to subset
    remove <- max_to_blast*2

    # if max_to_blast is more than the number of things to blast - take the number of things...
    if (length(input) < remove){
      remove <- length(input)
    }

    # make and write subset to blast file
    to_blast = input[(0:remove)]
    writeLines(to_blast, to_blast_path)

    # blast the subset
    output_table <- run_primer_blastn(to_blast_path, blast_db_path, ...)

    # add to previous results - if they exist
    append_table <- rbind(append_table, output_table)

    # remove duplicates
    append_table <- append_table %>% dplyr::group_by(saccver, sstart) %>% dplyr::filter(mismatch == min(mismatch))  %>% dplyr::distinct(saccver, sstart, .keep_all = TRUE)

    # update the file that contains primers to be blasted
    input = input[-(0:remove)]
    na.omit(input)
    writeLines(input, left_to_blast_path)

    #save results for later
    write.csv(append_table,
              file = append_table_path,
              row.names = FALSE)
  }

  unlink(left_to_blast_path)
  unlink(to_blast_path)

  append_table <- read.csv(append_table_path, colClasses = "character")

  # if output table is empty and give warning and stop
  if (nrow(append_table) <= 1){
    stop("No blast output generated.  Either no hits were found, or your compute environment could not support memory needs of the blastn step.  Try modifying parameters to reduce blast returns (e.g. align, max_to_blast, evalue, etc.)")

  }

  # parse amplicons from hits to forward and reverse hits
  # First isolate forward and reverse reads and rename columns
  F_only <-  dplyr::rename(dplyr::filter(append_table, grepl('forward',qseqid)), gi=sgi, accession = saccver, mismatch_forward = mismatch, forward_start = sstart, forward_stop = send )
  R_only <-  dplyr::rename(dplyr::filter(append_table, grepl('reverse',qseqid)), gi=sgi, accession = saccver, mismatch_reverse = mismatch, reverse_start = sstart, reverse_stop = send )


  # keep only the accessions with forward and reverse primer hits, and add a column for product length
  f_and_r <- dplyr::inner_join(x=F_only, y=R_only, by=c("accession", "gi", "staxids"))
  f_and_r <- dplyr::mutate(f_and_r, product_length=0)

  # calculate product length if F and R primer pairs are in correct orientation to make amplicon
  f_and_r <- dplyr::mutate(f_and_r, product_length = dplyr::case_when((forward_start < reverse_start & forward_start < forward_stop & reverse_stop < reverse_start ) ~ (as.numeric(reverse_start) - as.numeric(forward_start)),
                                                                      (forward_start > reverse_start & forward_start > forward_stop & reverse_stop > reverse_start) ~ (as.numeric(forward_start) - as.numeric(reverse_start)),))

  # remove all F and R primer pairs that would not make an amplicon
  f_and_r <- dplyr::filter(f_and_r, !is.na(product_length))

  # if table is empty and give warning and stop
  if (nrow( f_and_r ) <= 1){
    stop("No plausible amplicons were found.  Try modifying parameters to increase blast returns (e.g. num_fprimers_to_blast, num_rprimers_to_blast, align, evalue, etc.)")

  }

  #save unfiltered seeds output
  save_output_as_csv(f_and_r, "_unfiltered_get_seeds_local_output", out,
                     metabarcode_name)

  # keep only hits with acceptable product length
  f_and_r <- dplyr::filter(f_and_r, dplyr::between(product_length, minimum_length, maximum_length))

  # keep hits with accaptable number of mismatches
  f_and_r <- dplyr::filter(f_and_r, mismatch_forward <= mismatch & mismatch_reverse <= mismatch)

  # if table is empty and give warning and stop
  if (nrow( f_and_r ) <= 1){
    stop("Filtering removed all plausible amplicons.  Try modifying parameters to allow more amplicons to pass filter (e.g. minimum_length, maximum_length, mismatch, etc.)")
  }

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

  unlink(append_table_path)

  #return if you're supposed to
  if (return_table) {
    return(taxonomized_table)
  }
  else {
    return(NULL)
  }
}

`%>%` <- magrittr::`%>%`
