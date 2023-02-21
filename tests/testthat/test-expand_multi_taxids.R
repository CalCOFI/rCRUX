test_that("expand_multi_taxids works", {
  
  blast_db_path <- file.path(system.file(package = 'rCRUX', 'mock-db/blastdb'), 'mock-db')
  
  # fake data, random ids and accession from mock-db .. will it work?
  test_df <-
    data.frame(
      BLAST_db_taxids = c('86961;86962', '1117319'),
      accession = c('AB021891.1', 'KY987560.1')
    )
  
  result <- expand_multi_taxids(output_table = test_df, max_to_blast = 2, blast_db_path = blast_db_path)
  
  expect_equal(nrow(result), 2)
  
  
})
