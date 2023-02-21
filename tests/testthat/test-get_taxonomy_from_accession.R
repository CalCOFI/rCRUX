test_that("get_taxonomy_from_accession works", {
  
  input <-
    data.frame(
      accession = c('AB021891.1', #Anguilla australis schmidti
                    'AB021889.1') #Anguilla australis australis
    )
  
  accession_taxa_sql_path <- system.file(package = 'rCRUX', 'mock-db/taxonomizr-ncbi-db-small.sql')
  
  # hide taxonomizr warnings "cannot remove file '...', reason 'Permission denied'
  suppressWarnings(
    result <- get_taxonomy_from_accession(input = input, 
                                          accession_taxa_sql_path = accession_taxa_sql_path)
  )
  
  expected <-
    data.frame(
      accession = c("AB021891.1", "AB021889.1"),
      taxid = c(86962L, 86961L),
      species = c("Anguilla australis", "Anguilla australis"),
      superkingdom = c("Eukaryota", "Eukaryota"),
      kingdom = c("Metazoa", "Metazoa"),
      phylum = c("Chordata", "Chordata"),
      subphylum = c("Craniata", "Craniata"),
      superclass = c("Actinopterygii", "Actinopterygii"),
      class = c("Actinopteri", "Actinopteri"),
      subclass = c("Neopterygii", "Neopterygii"),
      order = c("Anguilliformes", "Anguilliformes"),
      family = c("Anguillidae", "Anguillidae"),
      subfamily = as.character(c(NA, NA)),
      genus = c("Anguilla", "Anguilla"),
      infraorder = as.character(c(NA, NA)),
      subcohort = as.character(c(NA, NA)),
      superorder = as.character(c(NA, NA)),
      superfamily = as.character(c(NA, NA)),
      tribe = as.character(c(NA, NA)),
      subspecies  = c("Anguilla australis schmidti" , "Anguilla australis australis"),
      subgenus = as.character(c(NA, NA)),
      species.group = as.character(c(NA, NA)),
      parvorder = as.character(c(NA, NA)),
      varietas = as.character(c(NA, NA))
    )
  
  expect_identical(result, expected = expected)
  
})
