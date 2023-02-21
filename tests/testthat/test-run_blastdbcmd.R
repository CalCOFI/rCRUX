test_that("run_blastdbcmd works on mock-db", {
  
  query_row <-
    data.frame(
      accession = 'HM485589.1',
      forward_stop = '668',
      reverse_stop = '1484'
    )
  
  blast_db_path <- file.path(system.file(package = 'rCRUX', 'mock-db/blastdb'), 'mock-db')
  
  result <- run_blastdbcmd(query_row = query_row, db = blast_db_path)
  
  # return accession with expanded by 1 range
  expect_true(grepl('>HM485589.1:667-1485', result[1]))
  
})

test_that("run_blastdbcmd fails expectedly on mock-db", {
  
  query_row <-
    data.frame(
      accession = 'XXXX',
      forward_stop = '668',
      reverse_stop = '1484'
    )
  
  blast_db_path <- file.path(system.file(package = 'rCRUX', 'mock-db/blastdb'), 'mock-db')
  
  expect_no_condition({
    result <- run_blastdbcmd(query_row = query_row, db = blast_db_path)
  })
  
  expect_true(attr(result, 'status') == 1)
  
})

