#' Create a small taxonomizr database for rCRUX testing
#' 
#' The script follows the `taxonomizr` db setup but without downloading the largest
#' table - only the names and nodes table (maybe only this is required?) That large
#' table - accession number:taxid - is manually added using `mock-db-sequences.map`.
#' The final local SQL DB will be saved as `taxonomizr-ncbi-db-small.sql`.
#' 
#' This script is expected to be run only during rCRUX development. The outputs
#' will can be used by end users if needed starting with `system.file(package = 'rCRUX', 'mock-db')`
#' 

library(magrittr)
library(dplyr)
library(RSQLite)

taxonomizr_db <- "inst/mock-db/taxonomizr-ncbi-db-raw.sql"

# Produces ~ 400 mb file
#taxonomizr::prepareDatabase(taxonomizr_db, getAccessions=FALSE, tmpDir = 'inst/mock-db/')

# Create a smaller version with out own accessions
small_db <- RSQLite::dbConnect(RSQLite::SQLite(), "inst/mock-db/taxonomizr-ncbi-db-small.sql")

# - Create accessionTaxa table
# need column names accession and taxa
accessionTaxa <-read.table('inst/mock-db/mock-db-sequences-taxid.map', sep = '\t', header = FALSE)
names(accessionTaxa) <- c('accession', 'taxa')
RSQLite::dbWriteTable(small_db, name = 'accessionTaxa', value = accessionTaxa)

# - Query big table
accession_taxa_db <- RSQLite::dbConnect(RSQLite::SQLite(), taxonomizr_db)
RSQLite::dbListTables(accession_taxa_db)

nodes <- dplyr::tbl(accession_taxa_db, 'nodes')
names <- dplyr::tbl(accession_taxa_db, 'names')

# ids to query nodes table
taxa_ids <-
  accessionTaxa$taxa

# node table
# traverse the tree to get all the taxonomy levels for a given ID
nodes_small <-
  lapply(taxa_ids, function(taxa_id){
    
    message(taxa_id)
    
    i = 1
    
    while (i == 1 || node_captured$parent != 1){
      
      node_captured <-
        nodes %>% 
        dplyr::filter(id == taxa_id) %>% 
        dplyr::collect()
      
      if( i == 1) {
        out <- node_captured
      } else {
        out <-
          rbind.data.frame(out, node_captured)
      }
      
      taxa_id <- node_captured$parent
      i <- i + 1
    }
    
    out
  })

nodes_small <- do.call(rbind.data.frame, nodes_small)
nodes_small <- dplyr::distinct(nodes_small)

# ids to query names table
nodes_ids <- unique(nodes_small$id, nodes_small$parent)

names_small <-
  names %>% 
  dplyr::filter(id %in% nodes_ids) %>% 
  dplyr::collect()

# update database
RSQLite::dbWriteTable(small_db, name = 'nodes', value = nodes_small)
RSQLite::dbWriteTable(small_db, name = 'names', value = names_small)

RSQLite::dbDisconnect(accession_taxa_db)
RSQLite::dbDisconnect(small_db)

# remove big db
#file.remove('inst/mock-db/taxonomizr-ncbi-db-raw.sql')
