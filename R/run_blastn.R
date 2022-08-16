#' Runs blastn with the input fasta as a query
#'
#' Takes the input string, writes it to a temporary file,
#' and calls blastn with that file as the query. Warning:
#' if a file at the path specified by temp already exists,
#' it will be overwritten then deleted.
#'
#' @param fasta a fasta-formatted string
#' @param temp: a file path to write a temporary fasta to
#' @param ncbi_bin: if not null use it as the parent directory for blastn
#' @return a tibble representing the blastn results
#' @export
run_blastn <- function(fasta, db_dir, temp = NULL, ncbi_bin = NULL,
                    evalue = 1e-6, align = 50000, coverage = 50, perID = 70) {
    # This is a hacky workaround to deal with the fact
    # that blastn wants a file path as a query
    # Ideally, we would find a way (perhaps a process substitution?)
    # to pass the string directly to blastn
    # The motivation for writing the fasta here rather than elsewhere is that
    # 1) This way, the script can still generally be visualized functionally
    # 2) It allows for easy changes if we ever figure out an elegant way to do
    # the handoff
    if (is.null(temp)) {
        temp <- tempfile()
    }

    message("Generated a temporary fasta at ", temp)

    writeLines(fasta, con = temp)

    # Determine arguments
    cores <- parallel::detectCores()

    message("Calling blastn. This may take a long time.")

    # System call
    if (is.null(ncbi_bin)) {
        blastn_output <- system2(command = "blastn",
                                    args = c("-db", db_dir,
                                    "-query", temp,
                                    "-outfmt", paste("\"6", "saccver", "length",
                                        "pident", "qacc", "slen", "sstart",
                                        "send", "sseq", "evalue", "staxids\""),
                                    "-evalue", evalue,
                                    "-num_alignments", align,
                                    "-qcov_hsp_perc", coverage,
                                    "-perc_identity", perID,
                                    "-num_threads ", cores),
                                    wait = TRUE,
                                    stdout = TRUE)
    }
    # Is there a way to add or remove the "/" based on need?
    # Maybe just look for a "/" on the end?
    else {
        blastn <- paste(ncbi_bin, "blastn", "/")
        blastn_output <- system2(command = blastn,
                                    args = c("-db", db_dir,
                                    "-query", temp,
                                    "-outfmt", paste("\"6", "saccver", "length",
                                        "pident", "qacc", "slen", "sstart",
                                        "send", "sseq", "evalue", "staxids\""),
                                    "-evalue", evalue,
                                    "-num_alignments", align,
                                    "-qcov_hsp_perc", coverage,
                                    "-perc_identity", perID,
                                    "-num_threads ", cores),
                                    wait = TRUE,
                                    stdout = TRUE)
    }

    file.remove(temp)

    # Format output
    column_names <-  c("accession",
                        "amplicon_length",
                        "pident",
                        "query_accession",
                        "accession_sequence_length",
                        "amplicon_start",
                        "amplicon_stop",
                        "sequence",
                        "evalue",
                        "BLAST_db_taxids")
    # This was a little hard for me to wrap my head around at first
    # as_tibble creates a one-column tibble with "value" as its col name
    # Curious if there is a better way to do this
    output_table <- blastn_output %>%
                    tibble::as_tibble() %>%
                    tidyr::separate(col = value, into = column_names,
                                    sep = "\t")
    return(output_table)
}

`%>%` <- magrittr::`%>%`