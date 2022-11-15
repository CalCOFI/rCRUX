#' Attach taxonomy data to an input table
#'
#' @description
#' This function will append taxids and taxonomy to any dataframe with a column
#' named accession that contains NCBI accessions (technically accession versions).
#'
#' @details
#' It takes a dataframe and searches for column called accession.  If present,
#' the column data is passed to [taxonomizr::accessionToTaxa()] to find the
#' corresponding taxids. The taxids are passed to [taxonomizr::getTaxonomy()]
#' and the following taxonomic information is retrieved for each of the following
#' ranks: species, superkingdom, kingdom, phylum, subphylum, superclass, class,
#' subclass, order, family, subfamily, genus, infraorder, subcohort, superorder,
#' superfamily, tribe, subspecies, subgenus, species group, parvorder, varietas.
#' New columns for taxid and each rank are appended to the dataframe.
#'
#' @param input a data.frame
#' @param accession_taxa_sql_path the path to an accessionTaxa sql
#'
#' @return the data.frame with taxonomy data
#'
#' @export
#'


get_taxonomizr_from_accession <- function(input, accession_taxa_sql_path,
                                        organize = TRUE) {
    if (!"accession" %in% colnames(input)) {
        stop("No `accession` column in input.")
    }
    input_taxids <- taxonomizr::accessionToTaxa(input$accession,
                                            accession_taxa_sql_path)

    input_taxonomy <- taxonomizr::getTaxonomy(input_taxids, accession_taxa_sql_path,
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
            Hint: this may be caused by 0-row inputs.")
    }

    if (organize) {
        # Arrange by taxonomy
        output <- output %>%
            dplyr::arrange(species) %>%
            dplyr::arrange(genus) %>%
            dplyr::arrange(family)  %>%
            dplyr::arrange(order) %>%
            dplyr::arrange(class) %>%
            dplyr::arrange(phylum) %>%
            dplyr::arrange(superkingdom)
    }
    return(output)
}
