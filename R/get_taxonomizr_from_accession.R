#' Attach taxonomy data to an input table
#'
#' @param input a data.frame
#' @param accessionTaxa_path the path to an accessionTaxa sql
#' @return the data.frame with taxonomy data
#' @export
get_taxonomizr_from_accession <- function(input, accessionTaxa_path,
                                        organize = TRUE) {
    if (!"accession" %in% colnames(input)) {
        stop("No `accession` column in input.")
    }
    input_taxids <- taxonomizr::accessionToTaxa(input$accession,
                                            accessionTaxa_path)

    input_taxonomy <- taxonomizr::getTaxonomy(input_taxids, accessionTaxa_path,
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