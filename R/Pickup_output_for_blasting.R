#' function to pickup_output_for_blasting new rounds
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



# Pick up where it left off
# Most simply, this does the "same thing" 5 times. function

pickup_output_for_blasting <- function(save_dir, output_table,
              unsampled_indices, too_many_ns, blastdbcmd_failed, num_rounds,
              blast_seeds_m) {

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

}
