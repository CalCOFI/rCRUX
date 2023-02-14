#' Combines matching fasta and taxonomy files, calculates product length, and if need adds taxids
#'
#' @description
#' combine_fasta_and_taxonomy takes a fasta file and a corresponding taxonomy file, and turns it
#' into an output .csv that can be used in [rCRUX::derep_and_clean_db()] and `rCRUX::compare_db()`.
#'
#'
#' @details
#' Fasta files are converted to a dataframe with columns accession and sequence. The product length
#' for each seqeunce is calculated. Taxonomy is joined to the corresponding accession, and taxonomic
#' ranks are separated by superkingdom, phylum, class, order, family, and genus and species. Lastly
#' taxonomy is added based on NCBI accession, or species id if the accession is not available.
#'
#' An rCRUX formatted summary .csv with the suffix "rCRUX_formatted_summary.csv" is generated in the
#' output directory indicated by output_directory_path.
#'
#'
#' @param output_directory_path the path to the output directory
#' @param metabarcode_name used to name the the files.
#' @param fasta_path the path to the fasta file of interest
#' @param taxonomy_path the path to the taxonomy file of interest
#' @param accession_taxa_sql_path the path to sql created by taxonomizr
#'        (e.g. accession_taxa_sql_path <- "/my/accessionTaxa.sql")
#' @param add.taxid does the user want to add taxid. The default is add.taxid = TRUE, which adds
#'        taxid bases on the status of NCBI.accession.
#' @param NCBI.accession is name of the sequence an the NCBI accession number.
#'        The default value is NCBI.accession = TRUE.  If NCBI.accession = FALSE
#'        the function will asign taxonomy using species name.
#'
#' @return NULL
#'
#' @export
#' @examples
#' 
#' \dontrun{
#'
#' output_directory_path <- "/my/directory/12S_mifish"
#' metabarcode_name <- "12S_mifish"
#' fasta_path = "/my/directory/12S_fasta_and_taxonomy/12S_.fasta"
#' taxonomy_path = "/my/directory/12S_fasta_and_taxonomy/12S_fasta_and_taxonomy/12S_taxonomy.txt"
#'
#' combine_fasta_and_taxonomy(
#'  output_directory_path = output_directory_path
#'  metabarcode_name = metabarcode_name,
#'  fasta_path = fasta_path, 
#'  taxonomy_path = taxonomy_path, 
#'  accession_taxa_sql_path = taxonomy_path
#')
#'}
#'


combine_fasta_and_taxonomy <-
  function(output_directory_path,
           metabarcode_name,
           fasta_path,
           taxonomy_path,
           accession_taxa_sql_path,
           NCBI.accession = TRUE,
           add.taxid = TRUE) {
    
  # turn a fasta file and taxonomy file into a dataframe with tax id and amplicon length
  
  # convert the fasta to a df get amplicon length
  fasta.df = phylotools::read.fasta(fasta_path)
  fasta.df <- dplyr::rename(fasta.df, accession = .data$seq.name, sequence = .data$seq.text)
  fasta.df <- dplyr::mutate(fasta.df, amplicon_length = nchar(sequence))
  
  # split taxonomic path
  tax <- utils::read.delim(taxonomy_path, header=FALSE)
  tax <- dplyr::rename(tax, accession = 1 , tax_path = 2)
  tax <- tidyr::separate(tax, col = 'tax_path', into = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"), ";")
  
  # join sequence and taxonomy
  full.df <- dplyr::left_join(fasta.df,tax, by = "accession", keep=FALSE)
  
  # add taxid from accession if ncbi accession is true, and from species name if accession is false - from species is a little less precise.
  if (NCBI.accession == TRUE & add.taxid == TRUE) {
    input_taxids <- taxonomizr::accessionToTaxa(full.df$accession, sqlFile = accession_taxa_sql_path)
    full.df <- dplyr::mutate(full.df, taxid = input_taxids)
    
  } else if (NCBI.accession == FALSE & add.taxid == TRUE) {
    taxid <- taxonomizr::getId(full.df$species,accession_taxa_sql_path)
    full.df <- dplyr::mutate(full.df, taxid_from_species = taxid)
    
  }
  
  utils::write.csv(full.df, 
                     file = file.path(output_directory_path, paste0(metabarcode_name, "_rCRUX_formatted_summary.csv" )), 
                     row.names = FALSE, col.names=TRUE)
  
}
