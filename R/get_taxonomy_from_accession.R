#' Attach taxonomy data to an input table using accession id
#'
#' @description
#' This function will append taxids and taxonomy to any dataframe with a column
#' named accession that contains NCBI accessions (technically accession versions).
#'
#' @details
#' It takes a dataframe and searches for a column called accession.  If present,
#' the column data is passed to [taxonomizr::accessionToTaxa()] to find the
#' corresponding taxids. The taxids are passed to [taxonomizr::getTaxonomy()]
#' and the following taxonomic information is retrieved for each of the following
#' ranks: species, superkingdom, kingdom, phylum, subphylum, superclass, class,
#' subclass, order, family, subfamily, genus, infraorder, subcohort, superorder,
#' superfamily, tribe, subspecies, subgenus, species group, parvorder, varietas.
#' New columns for taxid and each rank are appended to the dataframe.
#'
#' @param input a data.frame with column `accession`
#' @param accession_taxa_sql_path the path to local SQL database created by `taxonomizr`
#' @param arrange_taxonomy arrange rows of the output by taxonomy (superkingdom -> species)
#'
#' @return the data.frame with taxonomy data
#'
#' @export
#'
get_taxonomy_from_accession <- 
  function(input, 
           accession_taxa_sql_path,
           arrange_taxonomy = TRUE) {
    
    if (!file.exists(accession_taxa_sql_path)) {
      stop("accession_taxa_sql_path does not exist.\n",
           "The path to the taxonomizr SQL file cannot be found. ",
           "Please revise the path provided:\n", accession_taxa_sql_path)
    }
    
    if (!"accession" %in% colnames(input)) {
      stop("No `accession` column in input.")
    }
    
    input_taxids <- 
      taxonomizr::accessionToTaxa(input$accession,
                                  accession_taxa_sql_path)
    
    input_taxonomy <- 
      taxonomizr::getTaxonomy(input_taxids, accession_taxa_sql_path,
                              desiredTaxa = c("species", "superkingdom",
                                              "kingdom", "phylum", "subphylum", "superclass",
                                              "class", "subclass", "order", "family",
                                              "subfamily", "genus", "infraorder", "subcohort",
                                              "superorder", "superfamily", "tribe",
                                              "subspecies", "subgenus", "species group",
                                              "parvorder", "varietas"))
    
    output <-
      dplyr::mutate(input, taxid = input_taxids, data.frame(input_taxonomy))
    
    if (!"species" %in% colnames(output)) {
      stop("Failed to create column `species` in output.
            Hint: Is your data frame empty?")
    }
    
    if (arrange_taxonomy) {
      output <- 
        output %>%
        dplyr::arrange('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
    }
    
    output
  }

get_taxonomizr_from_accession <- function(){
  .Deprecated('get_taxonomy_from_accession', old = 'get_taxonomizr_from_accession')
}