# This line is to ensure this is testing the latest version of RCRUX.dev
devtools::install_github("LunaGal/RCRUX.dev")

accession_taxa_path <- "/data/home/galoscarleo/taxonomy/accessionTaxa.sql"
db_path <- "/data/home/galoscarleo/nt"


sink("short_test_out.txt")
err <- file("short_test_err.txt", open = "w")
sink(err, type = "message")
tryCatch(
    RCRUX.dev::get_blast_seeds("TAGAACAGGCTCCTCTAG", "TTAGATACCCCACTATGC",
        "short_test/", "12S_V5F1", accession_taxa_path,
        organism = c("7776", "7777"), return_table = FALSE),
    error = function(e) {
        message(e, "\n")
    }
)
closeAllConnections()

sink("large_test_out.txt")
err <- file("large_test_err.txt", open = "w")
sink(err, type = "message")
tryCatch(
    RCRUX.dev::get_blast_seeds("TAGAACAGGCTCCTCTAG", "TTAGATACCCCACTATGC",
                         "large_test/", "12S_V5F1", accession_taxa_path,
                         num_permutations = 25, hitsize = "1000000",
                         evalue = "100000", word_size = "6",
                         MAX_TARGET_PER_TEMPLATE = "5000",
                         NUM_TARGETS_WITH_PRIMERS = "500000",
                         organism = c("7776"), return_table = FALSE),
    error = function(e) {
        message(e, "\n")
    }
)
closeAllConnections()

sink("many_ns_test_out.txt")
err <- file("many_ns_test_err.txt", open = "w")
sink(err, type = "message")
tryCatch(
    RCRUX.dev::get_blast_seeds("GGWACWGGWTGAACWGTWTAYCCYCC",
                                "TANACYTCnGGRTGNCCRAARAAYCA",
                                "many_ns_test/", "CO1_063022",
                                accession_taxa_path,
                                num_permutations = 20, hitsize = "1000000",
                                evalue = "100000", word_size = "6",
                                MAX_TARGET_PER_TEMPLATE = "5",
                                NUM_TARGETS_WITH_PRIMERS = "500000",
                                organism = c("33208"), return_table = FALSE),
    error = function(e) {
        message(e, "\n")
    }
)
closeAllConnections()

sink("rcrux_blast_short_test_out.txt")
err <- file("rcrux_blast_short_test_err.txt", open = "w")
sink(err, type = "message")
tryCatch(
    RCRUX.dev::rcrux_blast("short_test/12S_V5F1/12S_V5F1_primerTree_output_with_taxonomy.csv",
                            "blast_test_save", db_path, accession_taxa_path),
    error = function(e) {
        message(e, "\n")
    }
)
closeAllConnections()

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