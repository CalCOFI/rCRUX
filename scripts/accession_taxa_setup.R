#Uses the taxonomizr library to generate an sqlite database

input <- blast
accessionTaxa_path <- 'taxonomizr_data/accessionTaxa.sql'
input_taxid <- accessionToTaxa(input$accession, accessionTaxa_path)

input_taxonomy <- getTaxonomy(input_taxid,accessionTaxa_path,desiredTaxa = c("species","superkingdom", "kingdom", "phylum", "subphylum", "superclass", "class", "subclass", "order", "family", "subfamily", "genus", "infraorder", "subcohort", "superorder", "superfamily", "tribe", "subspecies", "subgenus", "species group", "parvorder", "varietas"))

input_taxonomy <- cbind('accession'=input$accession, 'taxID'=input_taxid, input_taxonomy)
input_taxonomy <- as_tibble(input_taxonomy)
# Join the blast output and taxonomy tibbles
output <- full_join(input, input_taxonomy, by = "accession")