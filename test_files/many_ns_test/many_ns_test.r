accession_taxa_path <- "/data/home/galoscarleo/taxonomy/accessionTaxa.sql"
db_path <- "/data/home/galoscarleo/nt"

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