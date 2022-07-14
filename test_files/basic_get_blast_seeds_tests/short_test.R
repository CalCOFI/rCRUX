# Test the basic functionality of get_blast_seeds
sink("short_test.txt")
accession_taxa_path <- "/data/home/galoscarleo/taxonomy/accessionTaxa.sql"

RCRUX.dev::get_blast_seeds("TAGAACAGGCTCCTCTAG", "TTAGATACCCCACTATGC", "",
                "12S_V5F1", accession_taxa_path, organism = c("7776", "7777"),
                return_table = FALSE)
sink()