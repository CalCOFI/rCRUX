test_that("run_primer_blastn works", {
  
  # Create temporary fasta file with 16S primers
  temp_fasta <- tempfile(fileext = '.fasta')
  
  test_primers <-
    c('>primer_forward', 'AGAGGAGCGCGGAATTCC',
      '>primer_reverse', 'TACCTTGTTACGACTT')
  
  writeLines(test_primers, temp_fasta)
  
  blast_db_path <- file.path(system.file(package = 'rCRUX', 'mock-db/blastdb'), 'mock-db')
  
  # expect parameters passed on correctly 
  expect_message(regexp = '-num_alignments 100 ', 
                 object = { result <- run_primer_blastn(primer_fasta = temp_fasta, db = blast_db_path, align = 100) })
  
  expect_false(file.exists(temp_fasta))
  expect_true(inherits(result, 'data.frame'))
  expect_identical(names(result), c("qseqid", "sgi", "saccver", "mismatch", "sstart", "send", "staxids"))
  
})

test_that("run_primer_blastn returns correct errors", {
  
  temp_fasta <- tempfile(fileext = '.fasta')
  
  test_primers <-
    c('>primer_forward', 'AGAGGAGCGCGGAATTCC',
      '>primer_reverse', 'TACCTTGTTACGACTT')
  
  writeLines(test_primers, temp_fasta)
  
  #
  blast_db_path <- file.path(system.file(package = 'rCRUX', 'mock-db/blastdb'), 'mock-db123')
  # Looks for binaries first, if given
  expect_error(run_primer_blastn(primer_fasta = temp_fasta, db = blast_db_path, ncbi_bin = 'fail_me'), regexp = "The NCBI binaries could not be found")
  # Checks DB after finding binaries
  ncbi_bin <- 'C:/Program Files/NCBI/blast-BLAST_VERSION+/bin'
  if (dir.exists(ncbi_bin)){
    expect_error(run_primer_blastn(primer_fasta = temp_fasta, db = blast_db_path, ncbi_bin = ncbi_bin), regexp = "The BLAST DB could not be found")
  }
  expect_error(run_primer_blastn(primer_fasta = temp_fasta, db = blast_db_path), regexp = "The BLAST DB could not be found")
  
})

