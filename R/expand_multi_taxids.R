#' Expands returns from the nt database with multiple taxids

#' @description
#' The blast db downloaded from NCBIs FTP site has representative accessions
#' meaning identical sequences have been collapsed across multiple accessions
#' even if they have different taxid.
#' Here we identify representative accessions with multiple taxids, and unpack
#' all of the accessions that go into that representitive accessions.
#' Note - we are not identifying or unpacking the representative accessions
#' that report a single taxid
#'
#' The input table is similar to the summary.csv table produced by blast_seeds, but without taxonomy.
#' The output_table is in the same formate, but without multiple taxids in the blast taxids field
#'
#'
#' @param output_table table to use
#' @param max_to_blast passed to [rCRUX::blast_datatable()] and is the maximum
#'        number of entries to pass to blastcmdb.
#'        The default is max_to_blast = 1000 - the optimal number of reads to
#'        blast will depend on the user's environment (available RAM) and the
#'        number of possible hits (determined by marker and parameters)
#' @param blast_db_path a directory containing one or more blast-formatted database.
#'        For multiple blast databases, separate them with a space and add an extra set of quotes.
#'        (e.g blast_db_path <- "/my/ncbi_nt/nt" or blast_db_path <- '"/my/ncbi_nt/nt  /my/ncbi_ref_euk_rep_genomes/ref_euk_rep_genomes"')
#' @param ncbi_bin the directory that the blast+ suite is in. If NULL, the
#'        program will use your PATH environmental variable to locate them
#'
#' @return an output table
#'
#' @export
#' @examples
#' 
#' \dontrun{
#' expand_multi_taxids(output_table, max_to_blast)
#'}
expand_multi_taxids <- function(output_table, max_to_blast, blast_db_path, ncbi_bin = NULL){
  
  # Identify rows with multiple ids and filter into new dataframe
  multi_taxids <- dplyr::filter(output_table, grepl(';', .data$BLAST_db_taxids))
  
  # Exit function if no multi_taxids found
  if (nrow(multi_taxids) == 0){
    return(output_table)
  }
  
  # Continue if multi_taxids found ...
  message('Expanding multi taxids.')
  
  #remove multitaxids from output table
  clean_tax <- dplyr::setdiff(output_table, multi_taxids)
  
  #make a list of multi_taxid accessions for blastdbcmd
  accession <- paste0(multi_taxids$accession, collapse=",")
  
  
  # subset accessions and run through blastdbcmd.  This should solve memory problems
  list <- multi_taxids$accession
  expand_multi_taxids_output <- NULL
  
  # blastdbcmd system2 command value
  if (!is.null(ncbi_bin)){
    blastdbcmd <- file.path(ncbi_bin, 'blastdbcmd')
  } else {
    blastdbcmd = 'blastdbcmd'
  }
  
  while (length(list) > 0) {
    
    subset <- list[0:max_to_blast]
    subset<-subset[!is.na(subset)]
    find <- paste0(subset, collapse = ",")
    
    # send accessions through blastdbcmd
    
    expand_multi_taxids_out <- 
      system2(blastdbcmd,
              args = c("-db", blast_db_path,
                       "-dbtype", "nucl",
                       "-entry", find,
                       "-outfmt", '"%a %T"'),
              stdout = TRUE, stderr = FALSE)
    
    
    expand_multi_taxids_output <- c(expand_multi_taxids_out, expand_multi_taxids_output)
    list <- list[!(list %in% subset)]
  }
  
  #make df to store fixed multi taxid info
  accession_df <- dplyr::select(multi_taxids, accession)
  
  # make df for blastdbcmd output
  column_names <-  c("accession",
                     "BLAST_db_taxids")
  
  accession_output_table <- tibble::as_tibble(expand_multi_taxids_output)
  
  accession_output_table <-  
    tidyr::separate(accession_output_table, 
                    col = .data$value, 
                    into = column_names,
                    sep = " ")
  
  # add row_id for sorting later
  accession_output_table <- dplyr::mutate(accession_output_table, row_id = dplyr::row_number())
  
  # left join accessions that had mutli_taxids  with the blastcmd output - there wil be NA's so sort by row number to keep related blastdbcm output together - accessions that are identical
  accession_df <- dplyr::full_join(accession_df, accession_output_table,  by = "accession", keep = TRUE)
  accession_df <- dplyr::arrange(accession_df, .data$row_id)
  
  # fill NA's in columns with the initial accession value
  accession_df <- zoo::na.locf(zoo::na.locf(accession_df), fromLast=TRUE)
  
  # left join original multi_taxid table with the new expanded table -  fleshes out the info
  accession_df <- dplyr::left_join(accession_df, multi_taxids, by = c("accession.x" = "accession"))
  
  # remove unnecessary columns
  accession_df <- dplyr::select(accession_df, -c('accession.x', 'BLAST_db_taxids.y', 'row_id'))
  
  # change name of columns for concatonating with the clean data
  accession_df <- dplyr::rename(accession_df, accession = 'accession.y', BLAST_db_taxids = 'BLAST_db_taxids.x')
  
  # add the expanded blastdbcmd output with the single taxid table.
  output_table <- rbind(accession_df, clean_tax)
  
  return(output_table)
}
