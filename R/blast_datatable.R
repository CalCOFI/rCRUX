#' Wrapper and controller of the iterative blast search implemented by [rCRUX::run_blastdbcmd_blastn_and_aggregate_resuts]
#'
#' @description
#' #' Controls the iterative blast search implemented by
#'  [rCRUX::run_blastdbcmd_blastn_and_aggregate_resuts], cleans the output, and
#' adds taxonomy.
#' Given a datatable with the column names of the datatable returned by
#' [rCRUX::get_seeds_remote()], or [rCRUX::get_seeds_local()], uses a random
#' stratified sample based on taxonomic rank to iteratively process the data
#' The random sample entires are sent to
#' [rCRUX::run_blastdbcmd_blastn_and_aggregate_resuts()], which uses blastdbcmd
#' to convert entries into fasta files, passes them to blastn to query local
#' blast formatted databases with those sequences. It compiles the results
#' of blastn into a data.frame that it cleans and returns with taxonomy added
#' using [rCRUX::get_taxonomy_from_accession]. Additionally, it saves its
#' state as text files in a specified directory with each iteration, allowing
#' the user to restart an interrupted run of [rCRUX::blast_seeds()].
#'
#'
#' @details
#' blast_datatable uses [rCRUX::run_blastdbcmd_blastn_and_aggregate_resuts()] to
#' run blastdbcmd and blastn to find sequences. It randomly samples rows of the
#' seed table based on the taxonomic rank supplied by the user. The user can
#' specify how many sequences can be blasted simultaneously using max_to_blast.
#' The random sample will be subset for blasting. Once all of the seeds of the
#' random sample are processed, they are removed from the dataframe as are the
#' seeds found as blast hits. blast-datatable repeats this process or stratified
#' random sampling until there are fewer reads remaining than max_to_blast, at
#' which point it blasts all remaining seeds. The final aggregated results are
#' cleaned for multiple blast taxids, hyphens, and wildcards.
#'
#' Note:
#' The blast db downloaded from NCBIs FTP site has representative
#' accessions meaning identical sequences have been collapsed across multiple
#' accessions even if they have different taxid. Here we identify representative
#' accessions with multiple taxids, and unpack all of the accessions that were
#' collapsed into that representative accessions.
#'
#' We are not identifying or unpacking the representative accessions
#' that report a single taxid
#'
#' Saving data:
#' blast_datatable uses files generated in
#' [rCRUX::run_blastdbcmd_blastn_and_aggregate_resuts()] that store intermediate
#' results and metadata about the search to local files as it goes. This allows
#' the function to resume a partially completed blast, partially mitigating
#' the consequences of encountering an error or experiencing other interruptions.
#' Interruptions while blasting a subset of a random stratified sample will
#' result in a loss of the remaining reads of the subsample, and may decrease
#' overall blast returns. The local files are written to `save_dir` by
#' [rCRUX::save_state()]. Manually changing these files is not suggested as
#' it can change the behavior of blast_datatable.
#'
#' Restarting an interrupted [rCRUX::blast_seeds()] run:
#' To restart from an incomplete blast_datatable, submit the previous command
#' again. Do not modify the paths specified in the previous command, however
#' parameter arguments (e.g. rank, max_to_blast) can be modified. blast_datable
#' will automatically detect save files and resume from where it left off.
#'
#' Warning: If you are resuming from an interrupted blast, make sure you supply
#' the same data.frame for `blast_seeds`. If you intend to start a new blast,
#' make sure that there is not existing blast save data in the directory you
#' supply for `save_dir`.
#'
#' Note: blast_datatable does not save intermediate data
#' from blastdbcmd, so if it is interupted while getting building the fasta to
#' submit to blastn it will need to repeat some work when resumed. The argument
#' `max_to_blast` controls the frequency with which it calls blastn, so it can
#' be used to make blast_datatable save more frequently.
#'
#' @param blast_seeds a data.frame formatted like the output from
#'        get_seeds_remote or get_seeds_local
#' @param save_dir a directory in which to create files representing the
#'        current state
#' @param blast_db_path a directory containing one or more blast-formatted database.
#'        For multiple blast databases, separate them with a space and add an extra set of quotes.
#'        (e.g blast_db_path <- "/my/ncbi_nt/nt" or blast_db_path <- '"/my/ncbi_nt/nt  /my/ncbi_ref_euk_rep_genomes/ref_euk_rep_genomes"')
#' @param accession_taxa_sql_path a path to an sql created by
#'        [taxonomizr::prepareDatabase()]
#' @param ncbi_bin the directory that the blast+ suite is in. If NULL, the
#'        program will use your PATH environmental variable to locate them
#' @param force_db if true, try to use blast databases that don't appear to
#'        be blast databases
#' @param sample_size the number of entries to sample per rank
#'        before calling blastn - errors if not enough entries per rank (default = 1)
#' @param max_to_blast is the maximum number of entries to accumulate into a
#'        fasta before calling blastn (default = 1000)
#' @param wildcards a character vector representing the number of wildcards to
#'        discard (default = "NNNNNNNNNNNN")
#' @param rank the column representing the taxonomic rank to randomly sample (default = genus)
#' @param random_seed sets the random value generator for random stratified sampling. Change the default (random_seed = NULL) for reproducible results.
#' @param ... additional arguments passed to [rCRUX::run_blastdbcmd_blastn_and_aggregate_resuts()]
#' @inheritDotParams run_blastdbcmd_blastn_and_aggregate_resuts
#'
#' @return A data.frame representing the output of blastn
#'
#' @export


blast_datatable <- function(blast_seeds, save_dir, blast_db_path, accession_taxa_sql_path,
                            ncbi_bin = NULL, force_db = FALSE,
                            sample_size = 1, wildcards = "NNNNNNNNNNNN", rank = 'genus', max_to_blast = 1000, random_seed = NULL, ...) {

  check_blast_plus_installation(ncbi_bin = if('ncbi_bin' %in% names(list(...))) ncbi_bin else NULL)
  check_blast_db(blast_db_path)
  if (!file.exists(accession_taxa_sql_path)) {
    stop("accession_taxa_sql_path does not exist.\n",
         "The path to the taxonomizr SQL file cannot be found. ",
         "Please revise the path provided:\n", accession_taxa_sql_path)
  }

  # Default values for tracker variables
  num_rounds <- 1
  too_many_ns <- NULL
  blastdbcmd_failed <- NULL
  output_table <- NULL
  blast_seeds_m <- blast_seeds
  blast_seeds_m$blast_status <- "not_done"
  unsampled_indices <- seq_along(blast_seeds_m$accession)

  if (!is.null(random_seed)){
    set.seed(random_seed)
  }

  while (length(unsampled_indices) > 0) {

    # Get previous information if they exist
    unsampled_indicies_file <- file.path(save_dir, "unsampled_indices.txt")

    if (file.exists(unsampled_indicies_file)) {

      message('Previous unsampled indices exist, continuing from there.')

      rounds_path <- file.path(save_dir, "num_rounds.txt")
      num_rounds <- as.numeric(readLines(con = rounds_path))

      ns_path <- file.path(save_dir, "too_many_ns.txt")
      too_many_ns <- as.numeric(readLines(con = ns_path))

      blastdbcmd_failed_path <- file.path(save_dir, "blastdbcmd_failed.txt")
      blastdbcmd_failed <- as.numeric(readLines(con = blastdbcmd_failed_path))

      unsampled_indices_path <- file.path(save_dir, "unsampled_indices.txt")
      unsampled_indices <- as.numeric(readLines(con = unsampled_indices_path))

      output_table_path <- file.path(save_dir, "output_table.txt")
      output_table <- utils::read.csv(output_table_path, colClasses = "character")

      blast_seeds_m_path <- file.path(save_dir, "blast_seeds_passed_filter.txt")
      blast_seeds_m <- utils::read.csv(blast_seeds_m_path, colClasses = "character")

    }

    # Information about state of blast
    message("BLAST round: ", num_rounds)
    message('  ', length(unsampled_indices), " indices left to process.\n")

    # Update status of blast seeds by labeling all reads not in the upsampled
    # indicies list as "done"
    blast_seeds_m$blast_status[-unsampled_indices] <- "done"

    # Collect indices to blast
    # TODO unclear this part
    # if unsampled indices are greater than the max to blast (default n = 1000),
    # the blast seed table will be randomly sampled by taxonomic ranks
    if (length(unsampled_indices) <= max_to_blast) {
      sample_indices <- unsampled_indices
    }

    else if (rank == 'all') {

      sample_indices <- unsampled_indices

    }
    else {
      # if more indices than the max_to_blast are present
      # randomly select entries (default is n=1) for each rank then turn the
      # accession numbers into a vector
      # set random.seed for reproducible results

      seeds_by_rank_indices <-
        blast_seeds_m %>%
        dplyr::group_by(!!!rlang::syms(rank)) %>%
        dplyr::slice_sample(n = sample_size) %>%
        dplyr::filter(.data$blast_status == 'not_done') %>%
        dplyr::pull(.data$accession)
      # search the original output blast_seeds for the indices (row numbers) to
      # be used as blast seeds and make vector or sample indices
      sample_indices <- which(blast_seeds_m$accession %in% seeds_by_rank_indices)
    }

    # clean up messages
    if (length(unsampled_indices) > max_to_blast) {
      message(
        rank, " has ", length(sample_indices), " unique occurrences in the blast seeds data table.\n",
        "These may be subset ...\n"
      )
    } else {
      message("The number of unsampled indices is less than or equal to the maximum number to be blasted.\n")
    }

    # update unsampled_indices by removing the sample_indices from the list
    unsampled_indices <-
      unsampled_indices[!(unsampled_indices %in% sample_indices)]

    # run blast command, blastn, and aggregate the results based on the the value
    # max_to_blast. If there are fewer indices for a rank than the max_to_blast
    # it will run. If not the number of indices to be blasted for a rank will be
    # broken into the max_to_blast value.

    while (length(sample_indices) > 0 ){

      # Get previous information if they exist
      unsampled_indicies_file <- file.path(save_dir, "unsampled_indices.txt")

      if (file.exists(unsampled_indicies_file)) {

        rounds_path <- file.path(save_dir, "num_rounds.txt")
        num_rounds <- as.numeric(readLines(con = rounds_path))

        ns_path <- file.path(save_dir, "too_many_ns.txt")
        too_many_ns <- as.numeric(readLines(con = ns_path))

        blastdbcmd_failed_path <- file.path(save_dir, "blastdbcmd_failed.txt")
        blastdbcmd_failed <- as.numeric(readLines(con = blastdbcmd_failed_path))

        unsampled_indices_path <- file.path(save_dir, "unsampled_indices.txt")
        unsampled_indices <- as.numeric(readLines(con = unsampled_indices_path))

        output_table_path <- file.path(save_dir, "output_table.txt")
        output_table <- utils::read.csv(output_table_path, colClasses = "character")

        blast_seeds_m_path <- file.path(save_dir, "blast_seeds_passed_filter.txt")
        blast_seeds_m <- utils::read.csv(blast_seeds_m_path, colClasses = "character")

      }

      # TODO - this section could use better layout to only run run_blastdbcmd_blastn_and_aggregate_resuts
      # once, with some object that was created with the if else if elses e.g. `sample_indices_in`
      #
      # sample_indices == unsampled; sample_indices_in <- unsampled
      # sample_indices <= max_to_blast; sample_indices_in <- sample_indices
      # else;  sample_indices_in <- head(sample_indices, max_to_blast)
      #
      # How to update sample indices?

      if (length(sample_indices) == length(unsampled_indices)) {

        message('tmp - length(sample_indices) == length(unsampled_indices)\n')

        run_blastdbcmd_blastn_and_aggregate_resuts(sample_indices = unsampled_indices,
                                                   save_dir = save_dir,
                                                   blast_seeds_m = blast_seeds_m,
                                                   ncbi_bin = ncbi_bin,
                                                   db = blast_db_path,
                                                   too_many_ns = too_many_ns,
                                                   blastdbcmd_failed = blastdbcmd_failed,
                                                   unsampled_indices = unsampled_indices,
                                                   output_table = output_table,
                                                   wildcards = wildcards,
                                                   num_rounds = num_rounds,
                                                   ...)

        unsampled_indices <- unsampled_indices[!(unsampled_indices)]

        break

      } else if (length(sample_indices) <= max_to_blast) {

        message('tmp - length(sample_indices) <= max_to_blast\n')

        run_blastdbcmd_blastn_and_aggregate_resuts(sample_indices = sample_indices,
                                                   save_dir = save_dir,
                                                   blast_seeds_m = blast_seeds_m,
                                                   ncbi_bin = ncbi_bin,
                                                   db = blast_db_path,
                                                   too_many_ns = too_many_ns,
                                                   blastdbcmd_failed = blastdbcmd_failed,
                                                   unsampled_indices = unsampled_indices,
                                                   output_table = output_table,
                                                   wildcards = wildcards,
                                                   num_rounds = num_rounds,
                                                   ...)

        sample_indices <- sample_indices[!(sample_indices)]

        break

      } else {

        message('tmp - Subsetting sample_indices\n')

        # take chunks of the sample indices that are equivalent to max_to_blast
        subset <- utils::head(sample_indices, max_to_blast)

        run_blastdbcmd_blastn_and_aggregate_resuts(sample_indices = subset,
                                                   save_dir = save_dir,
                                                   blast_seeds_m = blast_seeds_m,
                                                   ncbi_bin = ncbi_bin,
                                                   db = blast_db_path,
                                                   too_many_ns = too_many_ns,
                                                   blastdbcmd_failed = blastdbcmd_failed,
                                                   unsampled_indices = unsampled_indices,
                                                   output_table = output_table,
                                                   wildcards = wildcards,
                                                   num_rounds = num_rounds,
                                                   ...)

        # update sample indices
        sample_indices <- sample_indices[!(sample_indices %in% subset)]
      }

    }

    rm(output_table)
    rm(too_many_ns)
    rm(blastdbcmd_failed)
    rm(blast_seeds_m)

  }

  # Clean up final datatable by removing any hyphens, updating amplicon_length
  # and removing reads with too many Ns

  output_table_path <- file.path(save_dir, "output_table.txt")
  output_table <- utils::read.csv(output_table_path, colClasses = "character")

  # Remove hyphens and reads with multiple Ns in output and recount amplicon length.
  output_table <- dplyr::mutate(output_table, sequence = gsub("-", "", sequence))

  too_many_ns <- dplyr::filter(output_table, grepl(wildcards, sequence))

  output_table <- dplyr::setdiff(output_table, too_many_ns)

  output_table <- dplyr::mutate(output_table, amplicon_length = nchar(sequence))


  #### The blast db downloaded from NCBIs FTP site has representative accessions
  # meaning identical sequences have been collapsed across multiple accessions
  # even if they have different taxid.
  # Here we identify representative accessions with multiple taxids, and unpack
  # all of the accessions that go into that representative accessions.
  # Note - we are not identifying or unpacking the representative accessions
  # that report a single taxid

  output_table <- expand_multi_taxids(output_table, max_to_blast, blast_db_path, ncbi_bin = NULL)

  output_table_taxonomy <-
    suppressWarnings(
      get_taxonomy_from_accession(output_table, accession_taxa_sql_path)
    )

  return(output_table_taxonomy)

}
