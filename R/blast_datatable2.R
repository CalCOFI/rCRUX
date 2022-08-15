# The idea here is to be a function like blast_datatable that takes any saved
# data as arguments for restarting.
# This makes the function less annoying and I can still wrap it in something
# else.
blast_datatable2 <- function(blast_seeds, save_dir, db_dir, accession_taxa_path,
                            sample_size = 1000, wildcards = "NNNN",
                            num_rounds = 0, too_many_ns = NULL,
                            not_in_db = NULL, output_table = NULL,
                            unsampled_indices = NULL) {
    # So this is basically gonna be the same as blast_datatable, right?
    # But we need special logic to initialise unsampled_indices
    # Nothing too fancy
    if (is.null(unsampled_indices)) {
        unsampled_indices <- seq_along(blast_seeds$accession)
    }
    # We can fill this in later
}