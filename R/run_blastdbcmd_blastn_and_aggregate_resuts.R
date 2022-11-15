#' runs blastdbcmd, blastn and aggregates the results
#'
#' @param sample_indices the indices to sample
#' @param save_dir a directory in which to create files representing the
#'        current state
#' @param blast_seeds_m blast seeds table but with blast status update
#' @param ncbi_bin passed to [rCRUX::run_blastdbcmd()] [rCRUX::run_blastn()] is
#'        the path to blast+ tools if not in the user's path.  Specify only if
#'        blastn and blastdbcmd  are not in your path.
#'        The default is ncbi_bin = NULL - if not specified in path do the
#'        following: ncbi_bin = "/my/local/blast+_folder".
#' @param too_many_ns a vector of indices that result
#'        in a fasta with too many Ns
#' @param db the type of blast db - e.g. nt
#' @param db_dir path to the blast db
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
#' @return NULL
#' @export



run_blastdbcmd_blastn_and_aggregate_resuts <- function(sample_indices = sample_indices,
        save_dir, blast_seeds_m, db, ncbi_bin = NULL, too_many_ns, db_dir, blastdbcmd_failed,
        unsampled_indices, output_table, wildcards, num_rounds, ...) {

        print("align =")
        print(align)


  # run blastdbcmd on each
  # sort results into appropriate buckets
  aggregate_fasta <- NULL
  message(paste("Running blastdbcmd on", length(sample_indices), "samples."))
  pb <- progress::progress_bar$new(total = length(sample_indices))



  for (index in sample_indices) {
    fasta <- suppressWarnings(run_blastdbcmd(blast_seeds_m[index, ], db, ncbi_bin))

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



    save_state(save_dir, output_table, unsampled_indices, too_many_ns,
               blastdbcmd_failed, num_rounds, blast_seeds_m)



    if (!is.character(aggregate_fasta)) {
      #message("aggregate_fasta has value ", aggregate_fasta)
      message("No useable accession numbers. Proceeding to next round.")

    }
    else {

      # run blastn and aggregate results
      blastn_output <- run_blastn(fasta=aggregate_fasta, db_dir=db, ncbi_bin=ncbi_bin)

      if(nrow(blastn_output) == 0 && length(unsampled_indices) > 0) {

      message(nrow(blastn_output), " blast hits returned.")

      stop(paste(" ", " ", "Either", " ",
      "1. blastn is having trouble blasting the number of seeds selected
              - reduce the value for max_to_blast", " ", "or", " ",
      "2. There were no hits returned for your blastn search because there were
    no valid matches in the database
             - try running the function again, and perhaps increase the
               max_to_blast."," ", sep="\n"))

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

      unsampled_indices <-
      unsampled_indices[!unsampled_indices %in% sample_indices]


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

   # add new blast round
   num_rounds <- num_rounds + 1

   # update files
   save_state(save_dir, output_table, unsampled_indices, too_many_ns,
                 blastdbcmd_failed, num_rounds, blast_seeds_m)


}
