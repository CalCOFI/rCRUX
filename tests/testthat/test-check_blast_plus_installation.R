test_that("check_blast_plus_installation throw error on fake binary path", {
  expect_error(check_blast_plus_installation(ncbi_bin = '1234'), regexp = 'The NCBI binaries could not be found at')
})
