#' function to dereplicate blast_seed results based on identical sequences and clean up reads with multiple taxids
#'
#' @param output_directory_path the path to the output directory
#' @param summary_path the path to the input file
#' @param metabarcode_name used to name the subdirectory and the files.
#' @return NULL
#' @export

# function to summarize the rcrux summary data by identical sequence and identify / collapse reads with multiple assignments for a given taxonomic rank.
# "We ain't too pretty, we ain't too proud" - Billy Joel, "Only the good die young"


derep_and_clean_db <- function(output_directory_path, summary_path, metabarcode_name) {

  out <- paste0(output_directory_path, "/derep_and_clean_db/")
  suppressWarnings(dir.create(out))


  # read rcrux blast summary.csv
  summary <- read.csv(summary_path)

  # get relevant df to work with
  summary <- dplyr::select(summary, accession, amplicon_length, sequence, taxid, superkingdom, phylum, class, order, family, genus, species)

  #remove hyphens from sequence
  summary <-  dplyr::mutate(summary, sequence = gsub("-", "", sequence))

  # save hits with no tax path
  no_path_summary <- dplyr::filter(summary, is.na(phylum) & is.na(class) & is.na(family) & is.na(genus))

  # remove hits with no tax path
  summary <- dplyr::filter(summary, !is.na(phylum) & !is.na(class) & !is.na(family) & !is.na(genus))


  # merge accessions and ranks for identical sequence
  phy_sum <- dplyr::summarize(dplyr::group_by(summary, sequence), accession = paste0(accession, collapse = ", "), amplicon_length = paste0(unique(amplicon_length), collapse = ", "), taxid = paste0(unique(taxid), collapse = ", "), superkingdom = paste0(unique(superkingdom), collapse = ", "), phylum = paste0(unique(phylum), collapse = ", "), class = paste0(unique(class), collapse = ", "), order = paste0(unique(order), collapse = ", "), family = paste0(unique(family), collapse = ", "), genus = paste0(unique(genus), collapse = ", "),   species = paste0(unique(species), collapse = ", "))


  #count number of accessions make new column
  phy_sum <- dplyr::mutate(phy_sum, num_of_accessions = (stringr::str_count(accession, ",") + 1 ))

  # remove , NA from ranks - they are most likely due to env seq or issues with sequence submission

    phy_sum <-  dplyr::mutate(phy_sum, superkingdom = gsub(", NA", "", superkingdom))
    phy_sum <-  dplyr::mutate(phy_sum, phylum = gsub(", NA", "", phylum))
    phy_sum <- dplyr::mutate(phy_sum, class = gsub(", NA", "", class))
    phy_sum <- dplyr::mutate(phy_sum, order = gsub(", NA", "", order))
    phy_sum <- dplyr::mutate(phy_sum, family = gsub(", NA", "", family))
    phy_sum <- dplyr::mutate(phy_sum, genus = gsub(", NA", "", genus))
    phy_sum <- dplyr::mutate(phy_sum, species = gsub(", NA", "", species))



  #Identify rows with multiple ids
  # subset dups:
  sub_dups <- dplyr::filter(phy_sum, grepl(', ', superkingdom) | grepl(', ', phylum) | grepl(', ', class) | grepl(', ', order) | grepl(', ', family) | grepl(', ', genus) | grepl(', ', species))


  #remove single taxonomy sequence
  clean_tax <- dplyr::setdiff(phy_sum, sub_dups)


  # change rank to NA if multiple names
  sub_dup_to_NA <- sub_dups

    sub_dup_to_NA <- dplyr::mutate(sub_dup_to_NA, superkingdom = gsub(".*, .*", "NA", superkingdom))
    sub_dup_to_NA <- dplyr::mutate(sub_dup_to_NA, phylum = gsub(".*, .*", "NA", phylum))
    sub_dup_to_NA <- dplyr::mutate(sub_dup_to_NA, class = gsub(".*, .*", "NA", class))
    sub_dup_to_NA <- dplyr::mutate(sub_dup_to_NA, order = gsub(".*, .*", "NA", order))
    sub_dup_to_NA <- dplyr::mutate(sub_dup_to_NA, family = gsub(".*, .*", "NA", family))
    sub_dup_to_NA <- dplyr::mutate(sub_dup_to_NA, genus = gsub(".*, .*", "NA", genus))
    sub_dup_to_NA <- dplyr::mutate(sub_dup_to_NA, species = gsub(".*, .*", "NA", species))


  write.csv(sub_dups,
            file = paste(out, "Sequences_with_multiple_taxonomic_paths.csv", sep = "/"),
            row.names = FALSE)

  write.csv(sub_dup_to_NA,
            file = paste(out, "Sequences_with_lowest_common_taxonomic_path_agreement.csv", sep = "/"),
            row.names = FALSE)


  write.csv(clean_tax,
            file = paste(out, "Sequences_with_single_taxonomic_path.csv", sep = "/"),
            row.names = FALSE)

  write.csv(no_path_summary,
            file = paste(out, "Sequences_with_mostly_NA_taxonomic_paths.csv", sep = "/"),
            row.names = FALSE)

  paths_to_summary_tables <- c(paste(out, "Sequences_with_lowest_common_taxonomic_path_agreement.csv", sep = "/"), paste(out, "Sequences_with_single_taxonomic_path.csv", sep = "/"))

  representative_fasta_and_taxonomy(paths_to_summary_tables, metabarcode_name, out)

}



# function to get the representative fasta and taxonomy files for derep_and_clean_db

representative_fasta_and_taxonomy <- function(paths_to_summary_tables, metabarcode_name, output_directory_path){

  # read in the file paths
  concat <- readr::read_csv(paths_to_summary_tables, show_col_types = FALSE)

  # grab the first instance of an accession to make representative
  concat <- dplyr::mutate(concat, rep_accession = purrr::map(strsplit(concat$accession, split = ","), 1) )

  # if only one accession give one fasta description and if more than one give another
  concat <- dplyr::mutate(concat, rep_accession_number = ifelse(concat$num_of_accessions == 1, paste0(concat$rep_accession), ifelse(concat$num_of_accessions > 1, paste0(concat$rep_accession,"_representative_of_", concat$num_of_accessions, "_identical_accessions"), 0)))


  # Write a fasta
  fasta <- character(nrow(concat) * 2)
  fasta[c(TRUE, FALSE)] <- paste0(">", concat$rep_accession_number)
  fasta[c(FALSE, TRUE)] <- concat$sequence
  writeLines(fasta, paste0(output_directory_path, "/", metabarcode_name, "_derep_and_clean.fasta"))


  # Taxonomy file format (tidyr and dplyr)
  taxa_table <-  dplyr::select(concat,rep_accession_number, superkingdom, phylum, class, order, family, genus, species)
  taxa_table <-tidyr::unite(taxa_table,taxonomic_path, superkingdom:species, sep = ";", remove = TRUE, na.rm = FALSE)
  taxa_table <-dplyr::slice(taxa_table,-1)

  # Write the thing
  taxa_table_path <- paste0(output_directory_path, "/", metabarcode_name, "_derep_and_clean_taxonomy.txt")
  write.table(taxa_table, file = taxa_table_path, row.names = FALSE, col.names=FALSE, sep = "\t")

  # Count distinct taxonomic ranks - includes NA
  tax_rank_sum <- dplyr::summarise_at(concat,c('superkingdom', 'phylum','class','order','family','genus','species'),dplyr::n_distinct)

  # Write output to blast_seeds_output
  tax_rank_sum_table_path <- paste0(output_directory_path, "/", metabarcode_name, "_derep_and_clean_unique_taxonomic_rank_counts.txt")
  write.table(tax_rank_sum, file = tax_rank_sum_table_path, row.names = FALSE, col.names=TRUE, sep = ",")

}
