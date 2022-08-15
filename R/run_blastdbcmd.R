#' Extracts query parameters from a row and calls blastdbcmd
#'
#' Given a row from a blast_seed formatted data.frame, extracts arguments
#' for blastdbcmd, then returns the output of the blastdbcmd call.
#'
#' @param query_row a row from get_blast_seeds
#' @param db a directory with a blast-formatted database
#' @param ncbi_bin if not null use it as the parent directory for blastn
#' @return a fasta-formatted character vector
#' @export
run_blastdbcmd <- function(query_row, db, ncbi_bin = NULL) {
    # Extract arguments
    accession <- query_row$accession
    forward <- query_row$forward_stop
    reverse <- query_row$reverse_stop

    # Massage forward and reverse
    if (forward < reverse) {
        forward <- forward + 1
        reverse <- reverse - 1
    }
    else {
        # Swap them
        temp <- forward
        forward <- reverse
        reverse <- temp

        # Tighten
        # This could be done in fewer lines but I expanded it for clarity
        forward <- forward + 1
        reverse <- reverse - 1
    }

    seq_range <- paste0(forward, "-", reverse)

    # System call
    if (is.null(ncbi_bin)) {
        fasta <- system2("blastdbcmd", args = c("-db", db,
                                                "-dbtype", "nucl",
                                                "-entry", accession,
                                                "-range", seq_range),
                                                stdout = TRUE, stderr = FALSE)
        return(fasta)
    }
    else {
        blastdbcmd <- paste(ncbi_bin, "blastdbcmd", sep = "/")
        fasta <- system2(blastdbcmd, args = c("-db", db,
                                                "-dbtype", "nucl",
                                                "-entry", accession,
                                                "-range", seq_range),
                                                stdout = TRUE, stderr = FALSE)
        return(fasta)
    }
}