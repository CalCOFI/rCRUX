#' Extracts seed amplicon sequence from blast databases using blastdbcmd and
#' primer start and end coordinates stored in the seeds table
#'
#' @description
#' Given a row from a blast_seed formatted data.frame, extracts arguments
#' for blastdbcmd, then returns the output of the blastdbcmd call.
#'
#' @param query_row a row from get_seeds_local or get_seeds_remote
#' @param db a directory with a blast-formatted database
#' @param ncbi_bin is the path to blast+ tools if not in the user's path.
#'        Specify only if blastn and blastdbcmd  are not in your path.
#'        The default is ncbi_bin = NULL - if not specified in path do the
#'        following: ncbi_bin = "/my/local/blast+_folder".
#' @return a fasta-formatted character vector
#' @export


run_blastdbcmd <- function(query_row, db, ncbi_bin = NULL) {
    # Extract arguments
    accession <- query_row$accession
    forward <- as.numeric(query_row$forward_stop)
    reverse <- as.numeric(query_row$reverse_stop)

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
        blastdbcmd_path <- paste(ncbi_bin, "blastdbcmd", sep = "/")
        print(blastdbcmd_path)
        fasta <- system2(command= blastdbcmd_path, args = c("-db", db,
                                                "-dbtype", "nucl",
                                                "-entry", accession,
                                                "-range", seq_range),
                                                stdout = TRUE, stderr = FALSE)
        return(fasta)
    }
}
