# A few examples of get_blast_seeds

# These file directories need to be changed to locations on your device

blast_seeds_parent <- "D:/blast_seeds_test"
accession_taxa_path <- "D:/taxonomizr_data/accessionTaxa.sql"

# This example tries to maximize the amount of information downloaded

test1 <- get_blast_seeds("TAGAACAGGCTCCTCTAG", "TTAGATACCCCACTATGC",
                         blast_seeds_parent, "12S_V5F1", accession_taxa_path,
                         num_permutations = 25, hitsize = "1000000",
                         evalue = "100000", word_size = "6",
                         MAX_TARGET_PER_TEMPLATE = "5000",
                         NUM_TARGETS_WITH_PRIMERS = "500000",
                         organism = c("7776"), return_table = FALSE)

# This example uses the default values to get a smaller response from NCBI

test2 <- get_blast_seeds("TAGAACAGGCTCCTCTAG", "TTAGATACCCCACTATGC",
                         blast_seeds_parent, "12S_V5F1", accession_taxa_path,
                         organism = c("7776", "7777"), return_table = FALSE)

# This demonstrates the error handling
# "dog" is actually a valid organism! Sadly (but usefully for demonstrating error handling), "dogpf" is not.

test3 <- get_blast_seeds("TAGAACAGGCTCCTCTAG", "TTAGATACCCCACTATGC",
                         blast_seeds_parent, "12S_V5F1", accession_taxa_path,
                         organism = c("7776", "7777", "dogpf"), return_table = FALSE)