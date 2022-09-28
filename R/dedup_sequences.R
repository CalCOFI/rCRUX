#' function to deduplicate results based on identical sequences
#'
#' @param output_dir the path to the output directory
#' @param input_file the path to the input file
#' @param metabarcode name of metabarcode

#' @return NULL
#' @export

#library(data.table)
#library(dplyr)

#metabarcode <- "12S_MiFish"
#output_dir <- "~/Dropbox/CRUX_2.0/12S_MiFish_3-10000/rcrux_blast_output/"
#input_file"~/Dropbox/CRUX_2.0/12S_MiFish_3-10000/rcrux_blast_output/summary.csv"

#summary <- read.csv(input_file)

# get smaller df to work with
#summary <- summary %>% select(accession, amplicon_length, sequence, phylum, class, order, family, genus, species)

#remove hyphens
#summary <-  dplyr::mutate(summary, sequence = gsub("-", "", sequence))

# ugly but effective
#phy_sum <- summary %>% group_by(sequence) %>% summarize(accession = paste0(accession, collapse = ", "), amplicon_length = paste0(unique(amplicon_length), collapse = ", "), phylum = paste0(unique(phylum), collapse = ", "), class = paste0(unique(class), collapse = ", "), order = paste0(unique(order), collapse = ", "), family = paste0(unique(family), collapse = ", "), genus = paste0(unique(genus), collapse = ", "),   species = paste0(unique(species), collapse = ", "))

#phy_sum_path <- paste0(output_dir, "/", metabarcode, "_phy_sum.txt")
#write.table(phy_sum, file = phy_sum_path, row.names = FALSE, col.names=FALSE, sep = "\t")
