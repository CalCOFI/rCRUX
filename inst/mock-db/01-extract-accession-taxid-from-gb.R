#' Get accession and tax ID from genebank entry
#' 
#' File was exported from NCBI, see notes

stopifnot('Expecting to be in the rCRUX root project folder' = basename(getwd()) == 'rCRUX')

library(xml2)

sequence.xml <- xml2::read_xml('inst/mock-db/sequence.fasta.xml')

accession_taxid <-
  data.frame(
    accession = sequence.xml |>
      xml2::xml_find_all('//TSeq_accver') |>
      xml2::xml_text(),
    taxid = sequence.xml |>
      xml2::xml_find_all('//TSeq_taxid') |>
      xml2::xml_text()
  )

accession_taxid |>
  write.table(file = 'inst/mock-db/mock-db-sequences-taxid.map',
              sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

