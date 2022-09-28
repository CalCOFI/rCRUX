#' Take a datatable and return the results of BLASTing it
#'
#' @description
#' Given a datatable with the column names of the datatable returned by
#' [RCRUX.dev::get_blast_seeds()], use blastdbcmd to convert entries into
#' fasta files, then uses blastn to query ncbi databases for those
#' sequences. It compiles the results of blastn into a data.frame that it
#' returns. Additionally, it saves its state as text files in a specified
#' directory with each iteration.
#'
#' @details
#' blast_datatable uses blastdbcmd and blastn to find sequences. It samples rows
#' from `blast_seeds` and uses blastdbcmd to find a seqence that corresponds to
#' the accession number and forward and reverse stops recorded in the table.
#' blastdbcmd outputs sequences as .fasta-formatted strings, which
#' blast_datatable concatenates into a multi-line fasta, then passes to blastn
#' as an argument. blast-datatable repeats this process until no rows remain,
#' aggregating the results in a single data.frame.
#'
#' # Saving data
#' blast_datatable writes intermediate results and metadata about the search to
#' local files as it goes. This allows the function to resume a partially
#' completed blast, mitigating the consequences of encountering an
#' error or experiencing other interruptions. The local files are written to
#' `save_dir` by [RCRUX.dev::save_state()]. Manually changing these files is not
#' suggested as it can change the behavior of blast_datatable. To start from an
#' incomplete blast_datatable, specify the same save_dir as the incomplete
#' blast. blast_datable will automatically detect save files and resume from
#' where it left off.
#'
#' Warning: If you are resuming from an interrupted blast, make sure you supply
#' the same data.frame for `blast_seeds`. If you intend to start a new blast,
#' make sure that there is not existing blast save data in the directory you
#' supply for `save_dir`.
#'
#' Note that blast_datatable does not save intermediate data
#' from blastdbcmd, so if it is interupted while getting building the fasta to
#' submit to blastn it will need to repeat some work when resumed. The argument
#' `sample_size` controls the frequency with which it calls blastn, so it can
#' be used to make blast_datatable save more frequently.
#'
#' @param blast_seeds a data.frame formatted like the output from
#'        get_blast_seeds
#' @param save_dir a directory in which to create files representing the
#'        current state
#' @param db a directory with a blast-formatted database
#' @param accession_taxa_path a path to an sql created by
#'        [taxonomizr::prepareDatabase()]
#' @param ncbi_bin the directory that the blast+ suite is in. If NULL, the
#'        program will use your PATH environmental variable to locate them
#' @param force_db if true, try to use blast databases that don't appear to
#'        be blast databases
#' @param sample_size the number of entries to sample per rank
#'        before calling blastn - errors if not enough entries per rank
#' @param max_to_blast is the maximum number of entries to accumulate into a
#'        fasta before calling blastn
#' @param wildcards a character vector representing the number of wildcards to
#'        discard
#' @param rank the column representing the taxonomic rank to sample
#' @return A data.frame representing the output of blastn
#' @export
blast_datatable <- function(blast_seeds, save_dir, db, accession_taxa_path,
                            ncbi_bin = NULL, force_db = FALSE,
                            sample_size = 1, wildcards = "NNNN", rank = 'genus', max_to_blast = 1000) {



  if (!(check_db(db) || force_db)) {
    stop(db, " is probably not a blast database.
         Use force_db = TRUE to try it anyway.")
  }

  # Default values for tracker variables
  num_rounds <- 1
  too_many_ns <- NULL
  blastdbcmd_failed <- NULL
  output_table <- NULL
  blast_seeds_m <- blast_seeds
  blast_seeds_m$blast_status <- "not_done"
  unsampled_indices <- seq_along(blast_seeds_m$accession)

######### testing
#  message("1")

  # Pick up where it left off
#  if (file.exists(paste(save_dir, "unsampled_indices.txt", sep = "/"))) {

#    rounds_path <- paste(save_dir, "num_rounds.txt", sep = "/")
#    num_rounds <- as.numeric(readLines(con = rounds_path))

#    ns_path <- paste(save_dir, "too_many_ns.txt", sep = "/")
#    too_many_ns <- as.numeric(readLines(con = ns_path))

#    blastdbcmd_failed_path <- paste(save_dir, "blastdbcmd_failed.txt", sep = "/")
#    blastdbcmd_failed <- as.numeric(readLines(con = blastdbcmd_failed_path))

#    unsampled_indices_path <-
#      paste(save_dir, "unsampled_indices.txt", sep = "/")

#    unsampled_indices <-
#      as.numeric(readLines(con = unsampled_indices_path))

#    output_table_path <- paste(save_dir, "output_table.txt", sep = "/")
#    output_table <- read.csv(output_table_path, colClasses = "character")


#    blast_seeds_m_path <- paste(save_dir, "blast_seeds_passed_filter.txt", sep = "/")
#    blast_seeds_m <- read.csv(blast_seeds_m_path, colClasses = "character")

#  }

  while (length(unsampled_indices) > 0) {


    ######### testing
      message("5")

      if (file.exists(paste(save_dir, "unsampled_indices.txt", sep = "/"))) {

        rounds_path <- paste(save_dir, "num_rounds.txt", sep = "/")
        num_rounds <- as.numeric(readLines(con = rounds_path))

        ns_path <- paste(save_dir, "too_many_ns.txt", sep = "/")
        too_many_ns <- as.numeric(readLines(con = ns_path))

        blastdbcmd_failed_path <- paste(save_dir, "blastdbcmd_failed.txt", sep = "/")
        blastdbcmd_failed <- as.numeric(readLines(con = blastdbcmd_failed_path))

        unsampled_indices_path <-
          paste(save_dir, "unsampled_indices.txt", sep = "/")
        unsampled_indices <-
          as.numeric(readLines(con = unsampled_indices_path))

        output_table_path <- paste(save_dir, "output_table.txt", sep = "/")
        output_table <- read.csv(output_table_path, colClasses = "character")


        blast_seeds_m_path <- paste(save_dir, "blast_seeds_passed_filter.txt", sep = "/")
        blast_seeds_m <- read.csv(blast_seeds_m_path, colClasses = "character")

      }


    # information about state of blast
    message(paste("BLAST round", num_rounds))
    message(paste(length(unsampled_indices), "indices left to process."))


    # update status of blast seeds by labeling all reads no in the upsampled
    # indicies list as "done"
    blast_seeds_m$blast_status[-unsampled_indices] <- "done"

    # collect indices to blast
    # if unsampled indices are greater than the max to blast (default n = 1000), the blast seed table will be randomly sampled by taxonomic ranks

    ######### testing
      message("3")

    if (length(unsampled_indices) <= max_to_blast) {
      sample_indices <- unsampled_indices
    }
    else  {

    ######### testing
      message("4")


      # if more indices than the max_to_blast are present
      # randomly select entries (default is n=1) for each rank then turn the
      # accession numbers into a vector

      seeds_by_rank_indices <- dplyr::pull(dplyr::filter(dplyr::slice_sample(dplyr::group_by(blast_seeds_m,!!!rlang::syms(rank)), n=sample_size), blast_status == 'not_done'), accession)

      # search the original output blast_seeds for the indices (row numbers) to
      # be used as blast seeds and make vector or sample indices
      sample_indices <- which(blast_seeds_m$accession %in% seeds_by_rank_indices)
    }

    message(paste(rank, "has", length(sample_indices), "unique occurrences in the blast seeds data table."))
    message(paste("These will be subsets into chunks of ", max_to_blast, " indices or fewer" ))

    # update unsampled_indices by removing the sample_indices from the list
    unsampled_indices <-
      unsampled_indices[!(unsampled_indices %in% sample_indices)]

    # save_state(save_dir, output_table, unsampled_indices, too_many_ns,
    #              blastdbcmd_failed, num_rounds, blast_seeds_m)


    # run blast command, blastn, and aggregate the results based on the the value
    # max_to_blast.  If there are fewer indices for a rank than the max_to_blast
    # it will run.  If not the number of indices to be blasted for a rank will be
    # broken into the max_to_blast value.

    while (length(sample_indices) > 0 ){

    ######### testing
      message("1.6.7")

      # Pick up where it left off
      if (file.exists(paste(save_dir, "unsampled_indices.txt", sep = "/"))) {

      rounds_path <- paste(save_dir, "num_rounds.txt", sep = "/")
      num_rounds <- as.numeric(readLines(con = rounds_path))

      ns_path <- paste(save_dir, "too_many_ns.txt", sep = "/")
      too_many_ns <- as.numeric(readLines(con = ns_path))

     blastdbcmd_failed_path <- paste(save_dir, "blastdbcmd_failed.txt", sep = "/")
      blastdbcmd_failed <- as.numeric(readLines(con = blastdbcmd_failed_path))

      unsampled_indices_path <-
          paste(save_dir, "unsampled_indices.txt", sep = "/")

      unsampled_indices <-
          as.numeric(readLines(con = unsampled_indices_path))

      output_table_path <- paste(save_dir, "output_table.txt", sep = "/")
      output_table <- read.csv(output_table_path, colClasses = "character")


      blast_seeds_m_path <- paste(save_dir, "blast_seeds_passed_filter.txt", sep = "/")
      blast_seeds_m <- read.csv(blast_seeds_m_path, colClasses = "character")

      }

      if (length(sample_indices) <= max_to_blast) {

      ######### testing
        message("6")

      run_blastdbcmd_blastn_and_aggregate_resuts(sample_indices, save_dir,
            blast_seeds_m, db, ncbi_bin = NULL, too_many_ns, db_dir,
            blastdbcmd_failed, unsampled_indices, output_table, wildcards,
            num_rounds)

      ######### testing
      message("6.5")
      save_state(save_dir, output_table, unsampled_indices, too_many_ns,
                          blastdbcmd_failed, num_rounds, blast_seeds_m)

      break

      }else{

      ######### testing
        message("7")

        # take chunks of the sample indices that are equivalent to max_to_blast
        subset <- head(sample_indices, max_to_blast)


        run_blastdbcmd_blastn_and_aggregate_resuts(subset, save_dir,
              blast_seeds_m, db, ncbi_bin = NULL, too_many_ns, db_dir,
              blastdbcmd_failed, unsampled_indices, output_table, wildcards,
              num_rounds)

        # update sample indices
        sample_indices <- sample_indices[!(sample_indices %in% subset)]
      }

    }

  ######### testing
    message("8")

    num_rounds <- num_rounds + 1

  ######### testing
    message("9")



  ######### testing
        message("4.5")
      #save_state(save_dir, output_table, unsampled_indices, too_many_ns,
      #            blastdbcmd_failed, num_rounds, blast_seeds_m)

  ######### testing
        message("4.7")

        rm(output_table)
        rm(unsampled_indices)
        rm(too_many_ns)
        rm(blastdbcmd_failed)
        rm(blast_seeds_m)

  }

  # If we get a taxid from blastn can we just use that?

  ######### testing
    message("10")

  output_table_taxonomy <-
    get_taxonomizr_from_accession(output_table, accession_taxa_path)
  return(output_table_taxonomy)
}

# True if the db is a blast database, false if it's not


check_db <- function(db, ncbi_bin = NULL) {
  if (is.null(ncbi_bin)) {
    try(system2("blastdbcmd", args = c("-db", db, "-info"), stdout = FALSE)) == 0
  } else {
    blastdbcmd <- paste0(ncbi_bin, "blastdbcmd")
    try(system2(blastdbcmd, args = c("-db", db, "-info"), stdout = FALSE)) == 0
  }
}
