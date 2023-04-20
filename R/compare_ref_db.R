#' compare_ref_db() is a Reference database comparison function
#'
#' @description
#'  Compare two reference databases, produces a Venn Diagram,
#'  Krona plot, and summary table of overlapping and mismatching accessions and associated taxonomy.
#'
#'
#'
#' Note:
#' Krona plot requires the installation of krona tools [krona tools](https://github.com/marbl/Krona/wiki/Installing) on your local computer.
#'
#'
#'
#' @param ref_db_1 a data.frame formatted like the output from
#'        blast_seeds() metabarcode_name_taxonomy.txt
#' @param ref_db_2 a data.frame formatted like the output from
#'        blast_seeds() metabarcode_name_taxonomy.txt
#' @param ref_db_1_name a text string of ref_db_1 name
#' @param ref_db_2_name a text string of ref_db_2 name
#' @param out_dir a directory in which to create summary files and figures
#' @param format_ref_db_1 indicates indicates if a user needs to format an output with separate
#'        columns for accession, superkingdom, phylum, class, order, family, genus,
#'        species, into a taxonomy.txt formatted file with two columns:
#'        accession and superkingdom;phylum;class;order;family;genus;species
#'        default value is FALSE (format_ref_db_2=FALSE)
#' @param format_ref_db_2 indicates if a user needs to format an output with separate
#'        columns for accession, superkingdom, phylum, class, order, family, genus,
#'        species, into a taxonomy.txt formatted file with two columns:
#'        accession and superkingdom;phylum;class;order;family;genus;species
#'        default value is FALSE (format_ref_db_2=FALSE).
#'
#' @return 1) A data.frame representing the overlapping accessions of both reference databases
#' 2) A data.frame representing the unique accessions in reference database 1
#' 3) A data.frame representing the unique accessions in reference database 2
#' 4) A venn diagram of accessions
#' 5) Krona plots of the species found only in reference database 1 and 2
#' @export
#'
#' @examples
#'
#' \dontrun{
#' ref_db_1_path <- "~/rCRUX_output/12S_parameters_test_1/get_seeds_local/12S_filtered_get_seeds_local_output_with_taxonomy.csv"
#' ref_db_2_path <- "~/rCRUX_output/12S_parameters_test_1/blast_seeds_output/12S_taxonomy.txt"
#' ref_db_1_name <- "get_seeds"
#' ref_db_2_name <- "blast_seeds"
#' out_dir <- "~/rCRUX_output/12S_parameters_test_1/comp_steps/"
#'
#' compare_ref_db(ref_db_1_path,
#'               ref_db_2_path,
#'               ref_db_1_name,
#'               ref_db_2_name,
#'               out_dir, format_ref_db_1=TRUE)
#'
#' ref_db_1_path <- "~/rCRUX_output/12S_parameters_test_1/blast_seeds_output/12S_taxonomy.txt"
#' ref_db_2_path <- "~/rCRUX_output/12S_parameters_test_2/blast_seeds_output/12S_taxonomy.txt"
#' ref_db_1_name <- "parameters_test_1"
#' ref_db_2_name <- "parameters_test_2"
#' out_dir <- "~/rCRUX_output/12S_parameters_tests/comp_parameters/"
#'
#' compare_ref_db(ref_db_1_path,
#'               ref_db_2_path,
#'               ref_db_1_name,
#'               ref_db_2_name,
#'               out_dir)
#' }
#'



compare_ref_db <- function(ref_db_1_path, ref_db_2_path,ref_db_1_name,ref_db_2_name,out_dir, format_ref_db_1=FALSE, format_ref_db_2=FALSE) {

  suppressWarnings(dir.create(out_dir))

# if input files are not in accession s;p;c;o;f;g;s format but have separate columns for accession s p c o f g s  make taxonomy formated files and save to out_dir

  if (isTRUE(format_ref_db_1)){
    output_to_taxonomy_file(ref_db_1_path, "ref_db_1", out_dir)
    ref_db_1_path <- file.path(out_dir, "ref_db_1_taxonomy.txt")
  }

  if (isTRUE(format_ref_db_2)){
    output_to_taxonomy_file(ref_db_2_path, "ref_db_2", out_dir)
    ref_db_2_path <- file.path(out_dir, "ref_db_2_taxonomy.txt")
  }


# import files.  A bit redundant if you just made them...

  ref_db_1 <- read.table(ref_db_1_path, header=F, sep = "\t") %>%  tibble::as_tibble() %>% dplyr::rename(accession=V1, sum.taxonomy=V2) %>%
    tidyr::separate(sum.taxonomy, into = c("superkingdom","phylum","class","order","family","genus","species"), sep=";")

  ref_db_2 <- read.table(ref_db_2_path, header=F, sep = "\t") %>%  tibble::as_tibble() %>% dplyr::rename(accession=V1, sum.taxonomy=V2) %>%
    tidyr::separate(sum.taxonomy, into = c("superkingdom","phylum","class","order","family","genus","species"), sep=";")

#remove duplicates if they exist

  ref_db_1 <- ref_db_1 %>% dplyr::distinct()

  ref_db_2 <- ref_db_2 %>% dplyr::distinct()




# identify taxonomic similarity and differences between datasets

  setdiff(ref_db_1$accession, ref_db_2$accession) -> in_1_not_2
  setdiff(ref_db_2$accession, ref_db_1$accession) -> in_2_not_1
  intersect(ref_db_2$accession, ref_db_1$accession) -> in_1_and_2
  union(ref_db_2$accession, ref_db_1$accession) -> in_1_or_2

#make summary table
  table_1 <- data.frame(length(in_1_not_2),
                        length(in_2_not_1),
                        length(in_1_and_2),
                        paste0(round(length(in_1_not_2)/length(in_1_or_2) *100,2),"%"),
                        paste0(round(length(in_2_not_1)/length(in_1_or_2) *100,2),"%"),
                        paste0(round(length(in_1_and_2)/length(in_1_or_2) *100,2),"%"))


  colnames(table_1) <- c(paste0(ref_db_1_name, " Unique Accessions"),paste0(ref_db_2_name," Unique Accessions"),"Overlapping Accessions",paste0("Perc Unique to ",ref_db_1_name),paste0("Perc Unique to ",ref_db_1_name),"Perc Overlap")
  write.table(table_1,file=paste0(out_dir,"Comparison_summary.txt"),quote=FALSE, sep="\t")

## Make Venn Diagram


  x <- list(ref_db_1_name=ref_db_1$accession,
            ref_db_2_name=ref_db_2$accession)

  names(x) <- c(ref_db_1_name,ref_db_2_name)

  p1 <- ggVennDiagram::ggVennDiagram(x) + ggplot2::scale_color_brewer(palette = "Paired") +
    ggplot2::theme(panel.background=ggplot2::element_rect(fill = "white",colour = "white"),
          plot.background = ggplot2::element_rect(fill = "white",colour = "white"),
          legend.key = ggplot2::element_rect(fill = "white"),
          legend.background = ggplot2::element_rect(fill = "white", colour="white"))
  p1
  ggplot2::ggsave(filename = paste0(out_dir,"Venn_diagram.png"),
         plot=p1,
         width = 10,
         height=6)


## Build All Method Phyloseq Object
  ### make pseudo ASV Reads
  ref_db_1 %>%
    dplyr::mutate(., db= ref_db_1_name, value=1) -> ref_db_1_db

  ref_db_2 %>%
    dplyr::mutate(., db= ref_db_2_name, value=1) -> ref_db_2_db

  dplyr::full_join(ref_db_1_db,ref_db_2_db) %>%
    tidyr::pivot_wider(names_from = db, values_from = value, values_fill = 0) -> combined_db

  combined_db[,-c((length(colnames(combined_db))-1):length(colnames(combined_db)))] %>%
    dplyr::distinct() -> hashes_unique


  hashes_unique <- hashes_unique %>% tibble::as_tibble(rownames = "number")
  hashes_unique <- hashes_unique %>% dplyr::mutate(number = paste0("taxon_",number))

  combined_db %>%
    dplyr::left_join(hashes_unique) %>%
    tidyr::pivot_longer(., cols = ref_db_1_name:ref_db_2_name, names_to="sample_Sample", values_to = "Detected")-> combined_data_long

  # add Metadata
  combined_data_long %>%
    dplyr::select(sample_Sample) %>% dplyr::distinct() %>% as.data.frame -> sampledata

  rownames(sampledata) <- sampledata$sample_Sample
  phyloseq::sample_data(sampledata) -> sampledata

  #add TAXonomy
  combined_data_long %>%
    dplyr::select(-sample_Sample, -Detected,-accession) %>%
    dplyr::select(superkingdom,phylum, class, order, family, genus, species,number) %>%
    dplyr::distinct() -> taxonomy_table


  taxonomy_table %>% as.matrix() -> taxonomy_table_mat
  rownames(taxonomy_table_mat) <- taxonomy_table$number

  TAX = phyloseq::tax_table(taxonomy_table_mat[,-length(colnames(taxonomy_table_mat))])

  # add ASV table
  combined_data_long %>%
    dplyr::select(number, sample_Sample, Detected) %>%
    tidyr::pivot_wider(names_from = "sample_Sample", values_from = "Detected") -> combined_wide

  combined_wide %>%
    dplyr::select(-number) %>% as.matrix() -> otu_table
  rownames(otu_table) <- combined_wide$number

  OTU = phyloseq::otu_table(otu_table, taxa_are_rows = TRUE)
  physeq_obj = phyloseq(OTU, TAX, sampledata)

# error if there are zero differences
  physeq_obj_diff = phyloseq::prune_taxa(taxa_sums(physeq_obj) < 2, physeq_obj)

  phyloseq::sample_data(physeq_obj_diff)$sample_Sample <- c(paste("Unique_to_", ref_db_1_name),paste("Unique_to_", ref_db_2_name))
  phyloseq::sample_names(physeq_obj_diff)<-c(paste("Unique_to_", ref_db_1_name),paste("Unique_to_", ref_db_2_name))

#plot Phyloseq object
  psadd::plot_krona(physeq_obj_diff,paste0(out_dir,"db_comparison_unique"),"sample_Sample",trim=T)
  i=1

  physeq_obj_int = phyloseq::prune_taxa(taxa_sums(physeq_obj) > 1, physeq_obj)
  names <- sample_names(physeq_obj_int)
  physeq_obj_int_2 = phyloseq::prune_samples(names[1], physeq_obj_int)

  phyloseq::sample_data(physeq_obj_int_2)$sample_Sample <- c("Overlapping_Accessions")
  phyloseq::sample_names(physeq_obj_int_2)<-c("Overlapping_Accessions")

  psadd::plot_krona(physeq_obj_int_2,paste0(out_dir,"db_comparison_overlap"),"sample_Sample",trim=T)

  return(table_1)

}


# function to turn a table with separate accession, s, p, o, c, f, g, s columns into rCRUX taxonomy file with two columns: accession and s;p;c;o;f;g;s

output_to_taxonomy_file <- function(table_path, metabarcode, out_dir ){

  output_table <- read.table(table_path, header=T, sep = ",")
  output_table <- output_table %>% dplyr::select('accession', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species') %>%
    tidyr::unite(col = 'taxonomic_path', 'superkingdom':'species', sep = ";", remove = TRUE, na.rm = FALSE) %>%
    dplyr::slice(-1)

  taxa_table_path <- file.path(out_dir, paste0(metabarcode, "_taxonomy.txt"))
  utils::write.table(output_table, file = taxa_table_path, row.names = FALSE, col.names=FALSE, sep = "\t")

}
