test_that("check_blast_db returns nothing with DB found", {
  blast_db_path <- file.path(system.file(package = 'rCRUX', 'mock-db/blastdb'), 'mock-db')
  expect_no_condition(check_blast_db(blast_db_path))
})

test_that("check_blast_db returns error when no DB found", {
  blast_db_path <- file.path(system.file(package = 'rCRUX', 'mock-db/blastdb'), 'mock-db123')
  expect_error(check_blast_db(blast_db_path), regexp = "The BLAST DB could not be found")
})
