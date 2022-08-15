#' Helper function to save the state of a running blast search
#'
#' @param save_dir the directory to save in
#' @param output_table the outputs generated so far
#' @param unsampled_indices a vector of indices not yet sampled
#' @param too_many_ns a vector of indices that result
#'        in a fasta with too many Ns
#' @param not_in_db a vector of indices not found in the local
#'        database
#' @param num_rounds the number of rounds so far
#' @return NULL
#' @export
save_state <- function(save_dir, output_table, unsampled_indices, too_many_ns,
                        not_in_db, num_rounds) {
    if (!dir.exists(save_dir)) {
        dir.create(save_dir)
    }
    write.csv(output_table,
            file = paste(save_dir, "output_table.txt", sep = "/"),
            row.names = FALSE)
    writeLines(as.character(unsampled_indices),
                con = paste(save_dir, "unsampled_indices.txt", sep = "/"))
    writeLines(as.character(too_many_ns),
                con = paste(save_dir, "too_many_ns.txt", sep = "/"))
    writeLines(as.character(not_in_db),
                con = paste(save_dir, "not_in_db.txt", sep = "/"))
    writeLines(as.character(num_rounds),
                con = paste(save_dir, "num_rounds.txt", sep = "/"))
}