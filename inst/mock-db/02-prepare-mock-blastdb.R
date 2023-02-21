#' Create a mock blast database for rCRUX testing
#' 
#' A few 16S, 18S and MiFish 12S locus sequences were manually selected from
#' NCBI and saved to `mock-db-sequences.fasta` and their accession number:taxid
#' saved to `mock-db-sequences.map` (two column, no column names, tab-delimited).
#' This script builds the mock db in the `inst/mock-db` directory of the package. 
#' It assumes working from project root directory, thus keep note of relative file 
#' paths.
#' 
#' This script is expected to be run only during rCRUX development. The outputs
#' will can be used by end users if needed starting with `system.file(package = 'rCRUX', 'mock-db')`
#' 

stopifnot('Expecting to be in the rCRUX root project folder' = basename(getwd()) == 'rCRUX')

#' Build DB
result <-
  system2('makeblastdb', 
          args =  c('-in',  'inst/mock-db/mock-db-sequences.fasta',
                    '-dbtype nucl',
                    '-parse_seqids',
                    '-title "Mock rCRUX source database"',
                    '-taxid_map', 'inst/mock-db/mock-db-sequences-taxid.map',
                    '-out inst/mock-db/blastdb/mock-db'
          )
  )

stopifnot('Previous system command returned an error result code' = result == 0)

# Test DB
# Nseqs
result_text <-
  system(
    'blastdbcmd -db inst/mock-db/blastdb/mock-db -info',
    intern = TRUE
  )

sequence_number <-
  strsplit(paste(result_text, collapse = ''), split = '\t')[[1]][2]

stopifnot('expecting 34 sequences in db' = grepl('^34', sequence_number))

# By id
result_text <-
  system(
    'blastdbcmd -db inst/mock-db/blastdb/mock-db -dbtype nucl -entry KY815345.1',
    intern = TRUE
  )

stopifnot('Returned value not does not have the entry value' = grepl('KY815345.1', paste(result_text, collapse = '')))

result_text <-
  system(
    'blastdbcmd -db inst/mock-db/blastdb/mock-db -dbtype nucl -entry X56576.1 -range 1-50',
    intern = TRUE
  )

stopifnot('Returned value not does not have the entry value' = grepl('X56576.1', paste(result_text, collapse = '')))

result_text <-
  system(
    'blastdbcmd -db inst/mock-db/blastdb/mock-db -dbtype nucl -entry LC091907.1 -range 1-50',
    intern = TRUE
  )

stopifnot('Returned value not does not have the entry value' = grepl('LC091907.1', paste(result_text, collapse = '')))




