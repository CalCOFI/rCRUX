#' Runs [rCRUX::run_blastdbcmd()], [rCRUX::run_blastn()], and aggregates and
#' saves the results
#'
#' @description
#' It uses [rCRUX::run_blastdbcmd()] to find a seed sequence that corresponds to
#' the accession number and forward and reverse stops recorded in the seeds table.
#' [rCRUX::run_blastdbcmd()] outputs sequences as .fasta-formatted strings, which
#' run_blastdbcmd_blastn_and_aggregate_resuts concatenates into a multi-line
#' fasta, then passes to [rCRUX::run_blastn()] as an argument. The output of
#' [rCRUX::run_blastn()] is de-replicated by accession, and only the longest
#' read per replicates is retained in the output table. The run state is saved
#' and passed back to [rCRUX::blast_datatable()].
#'
#'
#' @param sample_indices the indices to sample
#' @param save_dir a directory in which to create files representing the
#'        current state
#' @param blast_seeds_m blast seeds table but with blast status update
#' @param ncbi_bin passed to [rCRUX::run_blastdbcmd()] [rCRUX::run_blastn()] is
#'        the path to blast+ tools if not in the user's path.  Specify only if
#'        blastn and blastdbcmd  are not in your path.
#'        The default is ncbi_bin = NULL - if not specified in path do the
#'        following: ncbi_bin = "/my/local/ncbi-blast-2.10.1+/bin/".
#' @param too_many_ns a vector of indices that result
#'        in a fasta with too many Ns
#' @param db path to the blast db
#' @param blastdbcmd_failed a vector of indices not found in the local
#'        database
#' @param unsampled_indices the indices that need to be sampled
#' @param output_table the table of results
#' @param wildcards is a character vector that represents the minimum number
#'        of consecutive Ns the user will tolerate in a given seed or hit
#'        sequence. The default is wildcards = "NNNN"
#' @param num_rounds number of rounds of blast
#' @param blastdbcmd_failed the indicies not found in your blast db
#' @param ... additional arguments passed to [rCRUX::run_blastn()]
#' @inheritDotParams run_blastn
#'
#' @return NULL
#' @export

run_blastdbcmd_blastn_and_aggregate_resuts <-
  function(sample_indices,
           save_dir,
           blast_seeds_m,
           ncbi_bin = NULL,
           db,
           too_many_ns,
           blastdbcmd_failed,
           unsampled_indices,
           output_table,
           wildcards,
           num_rounds,
           ...) {

    # Run blastdbcmd on each sample index
    # sort results into appropriate buckets
    aggregate_fasta <- NULL
    message(paste("Running blastdbcmd on", length(sample_indices), "samples.\n"))
    pb <- progress::progress_bar$new(total = length(sample_indices))

    for (index in sample_indices) {

      fasta <-
        run_blastdbcmd(query_row = blast_seeds_m[index, ],
                       db = db,
                       ncbi_bin = ncbi_bin)

      # Check status (if attribute present) and if not 0 - blastdbcmd failed
      blastdbcmd_failed_status <-
        !is.null(attr(fasta, 'status')) && attr(fasta, 'status') != 0

      # No status returned if successful fasta returned, check wildcard
      has_too_many_ns <-
        is.null(attr(fasta, 'status')) && length(grep(wildcards, fasta)) > 0

      if (blastdbcmd_failed_status) {
        blastdbcmd_failed <- append(blastdbcmd_failed, index)
      }
      else if (has_too_many_ns) {
        too_many_ns <- append(too_many_ns, index)
      }
      else {
        aggregate_fasta <- append(aggregate_fasta, fasta)
      }
      pb$tick()
    }

    save_state(save_dir = save_dir,
               output_table = output_table,
               unsampled_indices = unsampled_indices,
               too_many_ns = too_many_ns,
               blastdbcmd_failed = blastdbcmd_failed,
               num_rounds = num_rounds,
               blast_seeds_m = blast_seeds_m)

    if (!is.character(aggregate_fasta)) {
      #message("aggregate_fasta has value ", aggregate_fasta)
      message("No useable accession numbers. Proceeding to next round.")
    }
    else {

      # run blastn and aggregate results
      blastn_output <-
        run_blastn(fasta = aggregate_fasta,
                   db = db,
                   ncbi_bin = ncbi_bin,
                   ...)

      if (nrow(blastn_output) == 0 && length(unsampled_indices) > 0) {

        stop(
          nrow(blastn_output), " blast hits returned.\n",
          "\nEither\n\n",
          "1. blastn is having trouble blasting the number of seeds selected\n",
          "   with the given set of parameters\n\n",
          "or\n\n",
          "2. There were no hits returned for your blastn search because there were no valid matches in the database\n",
          "  "
        )

      }
      else {

        message('  ', nrow(blastn_output), " blast hits returned.")

      }

      # remove accession numbers found by blast
      # this is not the most elegant way to do it but it's not the worst...
      in_output <- blast_seeds_m$accession %in% blastn_output$accession

      in_output_indices <- seq_along(blast_seeds_m$accession)[in_output]

      unsampled_indices <-
        unsampled_indices[!unsampled_indices %in% in_output_indices]

      unsampled_indices <-
        unsampled_indices[!unsampled_indices %in% sample_indices]

      blast_seeds_m$blast_status[-unsampled_indices] <- "done"


      # Add output to existing output
      if (is.null(output_table)) {

        output_table <- blastn_output

      }
      else {

        output_table <- tibble::add_row(output_table, blastn_output)

      }


      # Remove duplicated accessions, keeping the longest sequence
      output_table <-
        output_table %>%
        dplyr::group_by(.data$accession) %>%
        dplyr::filter(.data$amplicon_length == max(.data$amplicon_length)) %>%
        dplyr::filter(!(duplicated(.data$accession))) %>%
        dplyr::ungroup()

    }

    # report number of total unique blast hits
    message('  ', nrow(output_table), " unique blast hits after this round.\n")

    # add new blast round
    num_rounds <- num_rounds + 1

    # update files
    save_state(save_dir = save_dir,
               output_table = output_table,
               unsampled_indices = unsampled_indices,
               too_many_ns =too_many_ns,
               blastdbcmd_failed = blastdbcmd_failed,
               num_rounds = num_rounds,
               blast_seeds_m = blast_seeds_m)


  }
