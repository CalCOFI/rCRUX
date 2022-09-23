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

  # Pick up where it left off
  # This could be improved in a bunch of ways tbh
  # Most simply, this does the "same thing" 5 times. function
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

    output_table_path <- paste(save_dir, "output_table.txt", sep = "/")
    output_table <- read.csv(output_table_path, colClasses = "character")

    blast_seeds_m_path <- paste(save_dir, "blast_seeds_passed_filter.txt", sep = "/")
    blast_seeds_m <- read.csv(blast_seeds_m_path, colClasses = "character")

  }


  while (length(unsampled_indices) > 0) {
    #blast_seeds_m <- dplyr::filter(blast_seeds, !is.na(superkingdom) & !is.na(phylum) & !is.na(class) & !is.na(order))



    # information about state of blast
    message(paste("BLAST round", num_rounds))
    message(paste(length(unsampled_indices), "indices left to process."))

    # update status of blast seeds by labeling all reads no in the upsampled
    # indicies list as "done"
    blast_seeds_m$blast_status[-unsampled_indices] <- "done"

    # collect indices to blast
    # if unsampled indices are greater than the max to blast (default n = 1000), the blast seed table will be randomly sampled by taxonomic ranks

    if (length(unsampled_indices) <= max_to_blast) {
      sample_indices <- unsampled_indices
    }
    else  {

      # if more indices than the max_to_blast are present
      # randomly select entries (default is n=1) for each rank then turn the
      # accession numbers into a vector
      seeds_by_rank_indices <- dplyr::pull(dplyr::filter(dplyr::slice_sample(dplyr::group_by(blast_seeds_m,!!!rlang::syms(rank)), n=sample_size), blast_status == 'not_done'), accession)

      # search the original output blast_seeds for the indices (row numbers) to
      # be used as blast seeds and make vector or sample indices
      sample_indices <- which(blast_seeds_m$accession %in% seeds_by_rank_indices)
    }


    # update unsampled_indices by removing the sample_indices from the list
    unsampled_indices <-
      unsampled_indices[!(unsampled_indices %in% sample_indices)]

    # run blastdbcmd on each
    # sort results into appropriate buckets
    aggregate_fasta <- NULL
    message(paste("Running blastdbcmd on", length(sample_indices), "samples."))
    pb <- progress::progress_bar$new(total = length(sample_indices))

    for (index in sample_indices) {
      fasta <- run_blastdbcmd(blast_seeds_m[index, ], db, ncbi_bin)

      # Maybe in these cases we can just append directly to output?
      # room for improvement here...
      # Well, in cases where the command has a non-0 exit status,
      # system2 sometimes (always?) returns a character vector of length 0
      # This causes an error because there are no characters to check, so
      # the if has nothing to operate on. This kludgey `or` fixes that.
      if (length(fasta) == 0 || nchar(fasta) == 0) {
        blastdbcmd_failed <- append(blastdbcmd_failed, index)
      }
      else if (length(grep(wildcards, fasta)) > 0) {
        too_many_ns <- append(too_many_ns, index)
      }
      else {
        aggregate_fasta <- append(aggregate_fasta, fasta)
      }
      pb$tick()
    }

    if (!is.character(aggregate_fasta)) {
      #message("aggregate_fasta has value ", aggregate_fasta)
      message("No useable accession numbers. Proceeding to next round.")

    }

    else {
      # run blastn and aggregate results
      blastn_output <- run_blastn(fasta=aggregate_fasta, db_dir=db, ncbi_bin=ncbi_bin)

      if(nrow(blastn_output) == 0 && length(unsampled_indices) > 0) {
        message(nrow(blastn_output), " blast hits returned.")
        stop("Blastn having trouble blasting the number of seeds selected.  Try using a taxonomic rank with fewer unique groups, and reduce the value for max_to_blast")
      }
      else {
        message(nrow(blastn_output), " blast hits returned.")

      }

      # remove accession numbers found by blast
      # this is not the most elegant way to do it but it's not the worst...
      in_output <- blast_seeds_m$accession %in% blastn_output$accession
      in_output_indices <- seq_along(blast_seeds_m$accession)[in_output]

      unsampled_indices <-
        unsampled_indices[!unsampled_indices %in% in_output_indices]

      # Add output to existing output
      if (is.null(output_table)) {
        output_table <- blastn_output
      }
      else {
        output_table <- tibble::add_row(output_table, blastn_output)
      }

      # Remove duplicated accessions, keeping the longest sequence
      output_table <- output_table %>%
        dplyr::group_by(accession) %>%
        dplyr::filter(amplicon_length == max(amplicon_length)) %>%
        dplyr::filter(!(duplicated(accession)))
      output_table <- dplyr::ungroup(output_table)
    }
    # report number of total unique blast hits
    message(nrow(output_table), " unique blast hits after this round.")

    # save the state of the blast
    num_rounds <- num_rounds + 1
    save_state(save_dir, output_table, unsampled_indices, too_many_ns,
               blastdbcmd_failed, num_rounds, blast_seeds_m)
  }

  # If we get a taxid from blastn can we just use that?
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
