test_that("get_seeds_local returns stuff", {

  output_directory_path_top <- tempdir()
  output_directory_path <- file.path(output_directory_path_top, '16S')
  dir.create(output_directory_path, showWarnings = FALSE)
  
  # Source DBs
  blast_db_path <- file.path(system.file(package = 'rCRUX', 'mock-db/blastdb'), 'mock-db')
  accession_taxa_sql_path <- system.file(package = 'rCRUX', 'mock-db/taxonomizr-ncbi-db-small.sql')
  
  # Primers
  # Nitrospira F probe, Bacteria 1492R
  forward_primer_seq <- 'AGAGGAGCGCGGAATTCC'
  reverse_primer_seq <- 'TACCTTGTTACGACTT'
  
  metabarcode_name <- "Nitrospira"
  
  get_seeds_local(forward_primer_seq = forward_primer_seq,
                  reverse_primer_seq = reverse_primer_seq,
                  output_directory_path = output_directory_path,
                  metabarcode_name = metabarcode_name,
                  accession_taxa_sql_path = accession_taxa_sql_path,
                  blast_db_path = blast_db_path, 
                  minimum_length = 5, 
                  maximum_length = 900,
                  return_table = FALSE,
                  num_threads = 'max')
  
  # Tests
  get_seeds_local_files <- list.files(file.path(output_directory_path, 'get_seeds_local'), full.names = TRUE)
  expect_equal(length(get_seeds_local_files), 4)
})

test_that("ncbi_bin error works", {
  
  output_directory_path_top <- tempdir()
  output_directory_path <- file.path(output_directory_path_top, '16S')
  dir.create(output_directory_path, showWarnings = FALSE)
  
  # Source DBs
  blast_db_path <- file.path(system.file(package = 'rCRUX', 'mock-db/blastdb'), 'mock-db')
  accession_taxa_sql_path <- system.file(package = 'rCRUX', 'mock-db/taxonomizr-ncbi-db-small.sql')
  
  # Primers
  # Nitrospira F probe, Bacteria 1492R
  forward_primer_seq <- 'AGAGGAGCGCGGAATTCC'
  reverse_primer_seq <- 'TACCTTGTTACGACTT'
  
  metabarcode_name <- "Nitrospira"
  
  # A no existent ncbi_bin
  expect_error(
    get_seeds_local(forward_primer_seq = forward_primer_seq,
                    reverse_primer_seq = reverse_primer_seq,
                    output_directory_path = output_directory_path,
                    metabarcode_name = metabarcode_name,
                    accession_taxa_sql_path = accession_taxa_sql_path,
                    blast_db_path = blast_db_path, 
                    minimum_length = 5, 
                    maximum_length = 900,
                    return_table = FALSE,
                    num_threads = 'max',
                    ncbi_bin = 'NO_REAL_PATH'),
    regexp = 'The NCBI binaries could not be found at NO_REAL_PATH'
  )

})

test_that("ncbi_bin works", {
  
  output_directory_path_top <- tempdir()
  output_directory_path <- file.path(output_directory_path_top, '16S')
  dir.create(output_directory_path, showWarnings = FALSE)
  
  # Source DBs
  blast_db_path <- file.path(system.file(package = 'rCRUX', 'mock-db/blastdb'), 'mock-db')
  accession_taxa_sql_path <- system.file(package = 'rCRUX', 'mock-db/taxonomizr-ncbi-db-small.sql')
  
  # Primers
  # Nitrospira F probe, Bacteria 1492R
  forward_primer_seq <- 'AGAGGAGCGCGGAATTCC'
  reverse_primer_seq <- 'TACCTTGTTACGACTT'
  
  metabarcode_name <- "Nitrospira"
  
  # Extract NCBI path - assumes installed on runner for testing
  message(Sys.getenv('PATH'))
  path_items <- strsplit(Sys.getenv('PATH'), ';')[[1]]
  ncbi_path <- grep('blast', path_items, ignore.case = TRUE, value = TRUE)
  
  message('using [', ncbi_path, '] in testing')
  
  get_seeds_local(forward_primer_seq = forward_primer_seq,
                  reverse_primer_seq = reverse_primer_seq,
                  output_directory_path = output_directory_path,
                  metabarcode_name = metabarcode_name,
                  accession_taxa_sql_path = accession_taxa_sql_path,
                  blast_db_path = blast_db_path, 
                  minimum_length = 5, 
                  maximum_length = 900,
                  return_table = FALSE,
                  num_threads = 'max',
                  ncbi_bin = ncbi_path)
  
  # Tests
  get_seeds_local_files <- list.files(file.path(output_directory_path, 'get_seeds_local'), full.names = TRUE)
  expect_equal(length(get_seeds_local_files), 4)
  
})
