accession_taxa_path <- "D:/taxonomizr_data/accessionTaxa.sql"
blast_seeds_parent <- "D:/emily_test_output"
testCO1 <- get_blast_seeds("GGWACWGGWTGAACWGTWTAYCCYCC", "TANACYTCnGGRTGNCCRAARAAYCA",
                          blast_seeds_parent, "CO1_063022", accession_taxa_path,
                          num_permutations = 2, hitsize='1000000', evalue='100000',
                          word_size='6', MAX_TARGET_PER_TEMPLATE = '5',
                          NUM_TARGETS_WITH_PRIMERS ='500000',
                          organism = c("33208"), return_table = FALSE)