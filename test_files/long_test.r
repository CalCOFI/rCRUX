accession_taxa_path <- "/data/home/galoscarleo/taxonomy/accessionTaxa.sql"
db_path <- "/data/home/galoscarleo/nt"

RCRUX.dev::get_blast_seeds("TAGAACAGGCTCCTCTAG", "TTAGATACCCCACTATGC",
                            "large_test/", "12S_V5F1", accession_taxa_path,
                            num_permutations = 25, hitsize = "1000000",
                            evalue = "100000", word_size = "6",
                            MAX_TARGET_PER_TEMPLATE = "5000",
                            NUM_TARGETS_WITH_PRIMERS = "500000",
                            organism = c("7776"), return_table = FALSE)