# Hate the name
# The idea here is to return a vector of fasta-strings
# From blastdbcmd
generate_fastas <- function(X, database) {
    dbcmd <- function(accession, forward, reverse) {
        # So we're gonna call blastdbcmd much like we do in run_blastdbcmd
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

        # System call
        fasta <- system2("blastdbcmd", args = c("-db", database,
                                                "-dbtype", "nucl",
                                                "-entry", accession,
                                                "-range",
                                                paste0(forward, "-", reverse)),
                                                stdout = TRUE)
        return(fasta)
    }
    apply(X, 1, dbcmd)
}