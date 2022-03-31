#A few examples of get_blast_seeds

blast_seeds_parent <- "D:/blast_seeds"
accession_taxa_path <- "D:/accessionTaxa.sql"

test1 <- get_blast_seeds("TAGAACAGGCTCCTCTAG", "TTAGATACCCCACTATGC",
                blast_seeds_parent, "12S_V5F1_1", accession_taxa_path,
                organism = c("7776"), return_table = FALSE)

test2 <- get_blast_seeds("TAGAACAGGCTCCTCTAG", "TTAGATACCCCACTATGC",
                blast_seeds_parent, "12S_V5F1_2", accession_taxa_path,
                organism = c("7776", "7777"), return_table = FALSE)
