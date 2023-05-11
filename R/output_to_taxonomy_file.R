#' output_to_taxonomy_file
#' function to turn a table with separate accession, s, p, o, c, f, g, s columns 
#' into rCRUX taxonomy file with two columns: accession and s;p;c;o;f;g;s
#' 
#' @param table_path the path
#' @param metabarcode the code
#' @param out_dir the output path
#' 
#' @returns nothing
#' 
#' @export

output_to_taxonomy_file <- function(table_path, metabarcode, out_dir ){
  
  output_table <- utils::read.table(table_path, header=T, sep = ",")
  
  output_table <- output_table %>% dplyr::select('accession', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species') %>%
    tidyr::unite(col = 'taxonomic_path', 'superkingdom':'species', sep = ";", remove = TRUE, na.rm = FALSE) %>%
    dplyr::slice(-1)
  
  taxa_table_path <- file.path(out_dir, paste0(metabarcode, "_taxonomy.txt"))
  
  utils::write.table(output_table, file = taxa_table_path, row.names = FALSE, col.names=FALSE, sep = "\t")
  
  invisible(NULL)
  
}