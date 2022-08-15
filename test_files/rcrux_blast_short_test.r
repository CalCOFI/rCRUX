accession_taxa_path <- "/data/home/galoscarleo/taxonomy/accessionTaxa.sql"
db_path <- "/data/home/galoscarleo/nt/nt"
seeds_path <- "large_test/12S_V5F1/12S_V5F1_primerTree_output_with_taxonomy.csv"

RCRUX.dev::rcrux_blast(seeds_path, db_path, accession_taxa_path, "blast_test_save3", sample_size = 100) # nolint
