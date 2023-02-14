#' Helper function to save the state of a running blast search
#'
#' @details
#' Used by [rCRUX::run_blastdbcmd_blastn_and_aggregate_resuts()] to cache output
#' data after completed blast searches.
#'
#' @param save_dir the directory to save in
#' @param output_table the outputs generated so far
#' @param blast_seeds_m blast seeds table but with blast status update
#' @param unsampled_indices a vector of indices not yet sampled
#' @param too_many_ns a vector of indices that result
#'        in a fasta with too many Ns
#' @param blastdbcmd_failed a vector of indices not found in the local
#'        database
#' @param num_rounds the number of rounds so far
#' @return NULL
#' @export


save_state <- function(save_dir, output_table, unsampled_indices, too_many_ns,
                        blastdbcmd_failed, num_rounds, blast_seeds_m) {
    if (!dir.exists(save_dir)) {
        dir.create(save_dir, showWarnings = FALSE)
    }
    utils::write.csv(output_table,
            file = file.path(save_dir, "output_table.txt"),
            row.names = FALSE)
    utils::write.csv(blast_seeds_m,
            file = file.path(save_dir, "blast_seeds_passed_filter.txt"),
            row.names = FALSE)
    writeLines(as.character(unsampled_indices),
                con = file.path(save_dir, "unsampled_indices.txt"))
    writeLines(as.character(too_many_ns),
                con = file.path(save_dir, "too_many_ns.txt"))
    writeLines(as.character(blastdbcmd_failed),
                con = file.path(save_dir, "blastdbcmd_failed.txt"))
    writeLines(as.character(num_rounds),
                con = file.path(save_dir, "num_rounds.txt"))
}
