#' Get representative_fasta_and_taxonomy
#'
#' function to get the representative fasta and taxonomy files for derep_and_clean_db
#'
#' @param paths_to_summary_tables character vector with 2 values
#' @param metabarcode_name name
#' @param output_directory_path output path
#'
representative_fasta_and_taxonomy <- function(paths_to_summary_tables, metabarcode_name, output_directory_path){

  # read in the file paths
  concat <- readr::read_csv(paths_to_summary_tables, show_col_types = FALSE)

  # grab the first instance of an accession to make representative
  concat <- dplyr::mutate(concat, rep_accession = purrr::map(strsplit(concat$accession, split = ","), 1) )

  # if only one accession give one fasta description and if more than one give another
  concat <-
    dplyr::mutate(concat,
                  rep_accession_number = ifelse(.data$num_of_accessions == 1,
                                                paste0(.data$rep_accession),
                                                ifelse(.data$num_of_accessions > 1,
                                                       paste0(.data$rep_accession,"_representative_of_", .data$num_of_accessions, "_identical_accessions"), 0)
                                                )
                  )


  # Write a fasta
  fasta <- character(nrow(concat) * 2)
  fasta[c(TRUE, FALSE)] <- paste0(">", concat$rep_accession_number)
  fasta[c(FALSE, TRUE)] <- concat$sequence
  writeLines(fasta, file.path(output_directory_path, paste0(metabarcode_name, "_derep_and_clean.fasta")))

  # Taxonomy file format (tidyr and dplyr)
  taxa_table <-  dplyr::select(concat, 'rep_accession_number', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
  taxa_table <- tidyr::unite(taxa_table, col = 'taxonomic_path', 'superkingdom':'species', sep = ";", remove = TRUE, na.rm = FALSE)


  # Write the thing
  taxa_table_path <- file.path(output_directory_path, paste0(metabarcode_name, "_derep_and_clean_taxonomy.txt"))
  utils::write.table(taxa_table, file = taxa_table_path, row.names = FALSE, col.names=FALSE, sep = "\t")

  # Count distinct taxonomic ranks - includes NA
  tax_rank_sum <- dplyr::summarise_at(concat,c('superkingdom', 'phylum','class','order','family','genus','species'),dplyr::n_distinct)

  # Write output to blast_seeds_output
  tax_rank_sum_table_path <- file.path(output_directory_path, paste0(metabarcode_name, "_derep_and_clean_unique_taxonomic_rank_counts.txt"))
  utils::write.table(tax_rank_sum, file = tax_rank_sum_table_path, row.names = FALSE, col.names=TRUE, sep = "\t")


}
