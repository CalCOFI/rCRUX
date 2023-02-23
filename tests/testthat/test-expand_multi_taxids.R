test_that("expand_multi_taxids works", {
  
  blast_db_path <- file.path(system.file(package = 'rCRUX', 'mock-db/blastdb'), 'mock-db')
  
  test_df <-
    data.frame(
      BLAST_db_taxids = c('1794900;1798809;1857655;1857656;1857657', '1117319', '59291;1735697'),
      accession = c('KT851545.1', 'KY987560.1', 'JQ661395.1')
    )
  
  result <- expand_multi_taxids(output_table = test_df, max_to_blast = 2, blast_db_path = blast_db_path)
  
  expected_result <-
    data.frame(
      accession = c("KT851545.1", "KT851550.1", "KT851549.1", "KT851552.1", "KT851553.1", "JQ661395.1", "KT375339.1", "KY987560.1"),
      BLAST_db_taxids = c("1794900", "1798809", "1857657", "1857656", "1857655", "59291",   "1735697", "1117319")
    )
  
 expect_identical(result, expected = expected_result)
  
  
})
