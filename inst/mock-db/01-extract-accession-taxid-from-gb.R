#' Get accession and tax ID from genebank entry
#' 
#' File was exported from NCBI, see notes

stopifnot('Expecting to be in the rCRUX root project folder' = basename(getwd()) == 'rCRUX')

sequence.gb <- readLines('inst/mock-db/sequence.gb')

accession_taxid <-
  data.frame(
    version = sub(pattern = 'VERSION\\s{0,}',  
                  replacement = '' ,
                  x = sequence.gb[grep('^VERSION', sequence.gb)]),
    accession = sub(pattern = 'ACCESSION\\s{0,}([[:alnum:]]*).*',  
                    replacement = '\\1' ,
                    x = sequence.gb[grep('^ACCESSION', sequence.gb)]),
    taxid = gsub(pattern = ".*taxon:|\"", 
                 replacement = '',
                 x = sequence.gb[grep('taxon:', sequence.gb)])
  )

accession_taxid[c('version', 'taxid')] |>
  write.table(file = 'inst/mock-db/mock-db-sequences-taxid.map',
              sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
