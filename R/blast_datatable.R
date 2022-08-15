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
#' @param sample_size the number of entries to accumulate into a fasta before
#'        calling blastn
#' @param wildcards a character vector representing the number of wildcards to
#'        discard
#' @return A data.frame representing the output of blastn
#' @export
blast_datatable <- function(blast_seeds, save_dir, db, accession_taxa_path,
                            ncbi_bin = NULL,
                            sample_size = 1000, wildcards = "NNNN") {

    # Default values for tracker variables
    num_rounds <- 1
    too_many_ns <- NULL
    not_in_db <- NULL
    output_table <- NULL
    unsampled_indices <- seq_along(blast_seeds$accession)

    # Pick up where it left off
    # This could be improved in a bunch of ways tbh
    # Most simply, this does the "same thing" 5 times. function
    if (file.exists(paste(save_dir, "unsampled_indices.txt", sep = "/"))) {
        rounds_path <- paste(save_dir, "num_rounds.txt", sep = "/")
        num_rounds <- as.numeric(readLines(con = rounds_path))

        ns_path <- paste(save_dir, "too_many_ns.txt", sep = "/")
        too_many_ns <- as.numeric(readLines(con = ns_path))

        not_in_db_path <- paste(save_dir, "not_in_db.txt", sep = "/")
        not_in_db <- as.numeric(readLines(con = not_in_db_path))

        unsampled_indices_path <-
            paste(save_dir, "unsampled_indices.txt", sep = "/")
        unsampled_indices <-
            as.numeric(readLines(con = unsampled_indices_path))

        output_table_path <- paste(save_dir, "output_table.txt", sep = "/")
        output_table <- read.csv(output_table_path)
    }

    while (length(unsampled_indices) > 0) {
        # information about state of blast
        message(paste("BLAST round", num_rounds))
        message(paste(length(unsampled_indices), "indices left to process."))

        # sample some of them, removing them from the vector
        sample_indices <- smart_sample(unsampled_indices, sample_size)
        unsampled_indices <-
            unsampled_indices[!(unsampled_indices %in% sample_indices)]

        # run blastdbcmd on each
        # sort results into appropriate buckets
        aggregate_fasta <- NULL
        message(paste("Running blastdbcmd on", length(sample_indices), "samples."))
        pb <- progress::progress_bar$new(total = length(sample_indices))
        for (index in sample_indices) {
            fasta <- run_blastdbcmd(blast_seeds[index, ], db, ncbi_bin)
            # Maybe in these cases we can just append directly to output?

            # So this is somewhat atrocious. Why do we do it this way?
            # Well, in cases where the command has a non-0 exit status,
            # system2 sometimes (always?) returns a character vector of length 0
            # This causes an error because there are no characters to check, so
            # the if has nothing to operate on. This kludgey `or` fixes that.
            if (length(fasta) == 0 || nchar(fasta) == 0) {
                not_in_db <- append(not_in_db, index)
            }
            else if (length(grep(wildcards, fasta)) > 0) {
                too_many_ns <- append(too_many_ns, index)
            }
            else {
                aggregate_fasta <- append(aggregate_fasta, fasta)
            }
            pb$tick()
        }

        # run blastn and aggregate results
        blastn_output <- run_blastn(aggregate_fasta, db, ncbi_bin)

        # remove accesssion numbers found by blast
        # this is not the most elegant way to do it but it's not the worst...
        in_output <- blast_seeds$accession %in% blastn_output$accession
        in_output_indices <- seq_along(blast_seeds$accession)[in_output]
        # this message is to verify that I am doing this right
        message(length(in_output_indices),
                "indices were removed by the filtration step.")
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

        # save the state of the blast
        num_rounds <- num_rounds + 1
        save_state(save_dir, output_table, unsampled_indices, too_many_ns,
            not_in_db, num_rounds)
    }

    # If we get a taxid from blastn can we just use that?
    output_table_taxonomy <-
        get_taxonomizr_from_accession(output_table, accession_taxa_path)
    return(output_table_taxonomy)
}