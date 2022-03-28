#A few examples of get_blast_seeds

blast_seeds_parent <- "D:/blast_seeds"
accession_taxa_path <- "D:/accessionTaxa.sql"

get_blast_seeds("TAGAACAGGCTCCTCTAG", "TTAGATACCCCACTATGC",
                blast_seeds_parent, "12S_V5F1", accession_taxa_path,
                organism = c("7776"), MAX_TARGET_PER_TEMPLATE = 10, return_table = FALSE)
