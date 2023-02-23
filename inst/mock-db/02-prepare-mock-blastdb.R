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

# 312 unique seqs, 347 taxids
stopifnot('expecting 34 sequences in db' = grepl('^312', sequence_number))

# By id
result_text <-
  system(
    'blastdbcmd -db inst/mock-db/blastdb/mock-db -dbtype nucl -entry KY815345.1',
    intern = TRUE
  )

result_text

stopifnot('Returned value not does not have the entry value' = grepl('KY815345.1', paste(result_text, collapse = '')))

result_text <-
  system(
    'blastdbcmd -db inst/mock-db/blastdb/mock-db -dbtype nucl -entry KY815345.1 -outfmt "%a %T"',
    intern = TRUE
  )

result_text

stopifnot('Should be two hits' = length(result_text) == 2)

temp_fasta <- tempfile()

writeLines(con = temp_fasta, 
           text = c('>multitaxid', 
                    "GCCATCAGCTTTAACTAAACTTACACATGCGAGTATCCGCACCCCTGTGAGAATGCCCCGCTACCTCCTGATTGGAGACG
AGGAGCTGGCATCAGGCACAATCTTCTTTAGCCCACGACGCCTTGCTAAGCCACACCCCCAAGGGAATTCAGCAGTGATA
AATATTAAGCCATGAGTGAAAACTTGACTTAGTTAGTGTTTACAGGGCCGGTAAACCTCGTGCCAGCCACCGCGGTTAGA
CGAGAGACCCAAGTGGATGGCCCTCGGCGTAAAGAGTGGTTAAGATAGTCCCAAACTAAGGCCAAACGACTTTTTAGCTG
TTATACGCGCGAGAAGACATGAAGCCCAACTACGAAAGTGGCCTTAAGCCCCTGACCCCACGAAAGCTAGGATACAAACT
GGGATTAGATACCCCACTATGCCTAGCCCTAAACATTGATAACACCCTACCCAGTTTATCCGCCCGGGAACTACGAGCGT
CAGCTTGAAACCCAAAGGACTTGGCGGTGCTTTAGATCCACCTAGAGGAGCCTGTTCTAGAACCGATAATCCCCGTTCAA
CCTCACCTTTTCTTGCTTATTCCGCCTATATACCACCGTCGTCAGCTTACTCTTTGAGGAACTAGCCGTAAGCGCAACTG
GTACAACCTAAAACGCCAGGTCGAGGTGTAGCGAATGGAAAGGGAAGAAATGGGCTACATTCAGTAATAATAATGCATAC
GAACGATGTTCTGAAATAAACATCTGAAGGAGGATTTAGTAGTAAGTAGAGAGCAGAGTGCTCTACTGAAGCCGGCCCTG
AAGCGCGTACACACCGCCCGTCACTCTCCCCGAGCTAACCAAACACATTACTAATAAACAAAACTTGCAAAGGGGAGGC"
           ))

result_text <-
  system2(command = 'blastn',
          args = c("-db", "inst/mock-db/blastdb/mock-db",
                   "-query", temp_fasta,
                   "-outfmt", paste("\"6", "saccver", "length",
                                    "pident", "qacc", "slen", "sstart",
                                    "send", "evalue", "staxids\"")),
          wait = TRUE,
          stdout = TRUE)

result_text_staxids <-
  do.call(rbind.data.frame, 
          strsplit(result_text, split = '\t'))[,9]

stopifnot('Should be mutlitaxids with ; ' = any(grepl(';', result_text_staxids)))
