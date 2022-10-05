#' function to deduplicate results based on identical sequences
#'
#' @param output_dir the path to the output directory
#' @param summary_path the path to the input file

#' @return NULL
#' @export

# function to summarize the rcrux summary data by identical sequence and identify / collapse reads with multiple assignments for a given taxonomic rank.
# "We ain't too pretty, we ain't too proud" - Billy Joel, "Only the good die young"


dedup <- function(output_dir, summary_path, rank = c("phylum", "class", "order", "family", "genus", "species"), remove_NA = TRUE ) {

  # read rcrux blast summary.csv
  summary <- read.csv(summary_path)

  # get relevant df to work with
  summary <- dplyr::select(summary, accession, amplicon_length, sequence, phylum, class, order, family, genus, species)

  #remove hyphens from sequence
  summary <-  dplyr::mutate(summary, sequence = gsub("-", "", sequence))

#ok to here

  # merge accessions and ranks for identical sequence
  #phy_sum <- summary %>% group_by(sequence) %>% summarize(accession = paste0(accession, collapse = ", "), amplicon_length = paste0(unique(amplicon_length), collapse = ", "), phylum = paste0(unique(phylum), collapse = ", "), class = paste0(unique(class), collapse = ", "), order = paste0(unique(order), collapse = ", "), family = paste0(unique(family), collapse = ", "), genus = paste0(unique(genus), collapse = ", "),   species = paste0(unique(species), collapse = ", "))

  phy_sum <- dplyr::summarize(group_by(summary, sequence), accession = paste0(accession, collapse = ", "), amplicon_length = paste0(unique(amplicon_length), collapse = ", "), phylum = paste0(unique(phylum), collapse = ", "), class = paste0(unique(class), collapse = ", "), order = paste0(unique(order), collapse = ", "), family = paste0(unique(family), collapse = ", "), genus = paste0(unique(genus), collapse = ", "),   species = paste0(unique(species), collapse = ", "))


  #count number of accessions make new column
  phy_sum <- dplyr::mutate(phy_sum, num_of_accessions = (stringr::str_count(accession, ",") + 1 ))

  # remove , NA from ranks - they are most likely due to env seq or issues with sequence submission

  if (remove_NA == TRUE){
    for (r in rank){
      phy_sum <- dplyr::mutate(phy_sum, !! r := sub(", NA", "", !! sym(r)))
    }

  }


  #Identify rows with multiple ids
  # subset dups:
  sub_dups <- dplyr::filter(phy_sum, grepl(', ', phylum) | grepl(', ', class) | grepl(', ', order) | grepl(', ', family) | grepl(', ', genus) | grepl(', ', species))


  #remove single taxonomy sequence
  clean_tax <- dplyr::setdiff(phy_sum, sub_dups)


  # change rank to NA if multiple names
  sub_dup_to_NA <- sub_dups

  for (r in rank){
    sub_dup_to_NA <- dplyr::mutate(sub_dup_to_NA, !! r := sub(".*, .*", "NA", !! sym(r)))
  }


  # write output

  write.csv(sub_dups,
            file = paste(output_dir, "References_with_multiple_taxonomic_ranks.csv", sep = "/"),
            row.names = FALSE)

  write.csv(sub_dup_to_NA,
            file = paste(output_dir, "References_with_NA_instead_of_multiple_taxonomic_ranks.csv", sep = "/"),
            row.names = FALSE)


  write.csv(clean_tax,
            file = paste(output_dir, "References_with_unique_taxonomic_ranks.csv", sep = "/"),
            row.names = FALSE)

}
