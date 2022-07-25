accession_taxa_path <- "/data/home/galoscarleo/taxonomy/accessionTaxa.sql"
db_path <- "/data/home/galoscarleo/nt"

sink("rcrux_blast_large_test_out.txt")
err <- file("rcrux_blast_large_test_err.txt", open = "w")
sink(err, type = "message")
tryCatch(
    RCRUX.dev::rcrux_blast("large_test/12S_V5F1/12S_V5F1_primerTree_output_with_taxonomy.csv",
                            "blast_test_save", db_path, accession_taxa_path),
    error = function(e) {
        message(e, "\n")
    }
)
closeAllConnections()