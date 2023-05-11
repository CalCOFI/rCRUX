#' reformats the output of [rCRUX::blast_seeds()]
#' and writes it to disk in a format compatible with Qiime2
#'
#' @description
#'  Converts an [rCRUX::blast_seeds()] formatted database into a Qiime2
#'  formatted database.
#'
#'  The functions reads the file into R, reformats it into Qiime desired format
#'  and writes it in the desired path.
#'
#' @details
#'
#'
#'  The input file will have two columns separated by a tab, and no headers. The
#'  first column contains the accession of the sequences, the second column
#'  includes the taxonomical annotation of those sequences, with the annotations
#'  superkingom, phylum, class, order, family, genus and species in a single
#'  string separated by colons.
#'
#'  The function reads in the input file using [utils::read.table()], separates
#'  the two fields and follows by separating the taxonomical fields into
#'  columns.
#'
#'   With the option `format_ref_db_1` set to `TRUE` (the default), the function
#'  expects a comma-separated file with separate columns for accession and the
#'  following taxonomical ranks: superkingdom, phylum, class, order, family,
#'  genus, species. The input database is reformatted using
#'  [rCRUX::output_to_taxonomy_file].
#'
#'  One the input file is in the desired format, the function reformats the
#'  taxonomical fields by adding the prefixes that identified each taxonomical
#'  rank.
#'
#'
#'  The function returns no value into the R session, but writes to file on
#'  disk.
#'
#' @param ref_db_1_path The path to a tab-separated file, formatted like the
#'  output from [rCRUX::blast_seeds()]
#' @param ref_db_1_name a text string of ref_db_1 name
#' @param format_ref_db_1 logical. If TRUE, assumes the input file comes in the
#'        format of a comma-separated file with named columns for accession,
#'        superkingdom, phylum, class, order, family, genus, species.
#'        Defaults to FALSE (format_ref_db_1=FALSE)
#' @param out_dir a path to the directory where to write the output database
#'
#' @returns NULL
#'
#' @examples
#' \dontrun{
#' convert_to_qiime2(
#'   ref_db_1_path = here(
#'     "eecurd_updated_dbs_20230420",
#'     "expanded_12S_20230420",
#'     "12S",
#'     "blast_seeds_output",
#'     "12S_taxonomy.txt"
#'   ),
#'   ref_db_1_name = "12S",
#'   out_dir = here()
#' )
#' }
#'
#' @export

convert_to_qiime2 <-
  function(ref_db_1_path,
           ref_db_1_name,
           out_dir,
           format_ref_db_1 = FALSE) {
    
    dir.create(out_dir, showWarnings = FALSE)
    
    # if input files are not in accession s;p;c;o;f;g;s format but have separate columns for accession s p c o f g s  make taxonomy formatted files and save to out_dir
    
    if (isTRUE(format_ref_db_1)) {
      output_to_taxonomy_file(ref_db_1_path, "ref_db_1", out_dir)
      ref_db_1_path <- file.path(out_dir, "ref_db_1_taxonomy.txt")
    }
    
    # import files.  A bit redundant if you just made them...
    
    ref_db_1 <-
      utils::read.table(ref_db_1_path, header = F, sep = "\t") %>%
      tibble::as_tibble() %>% 
      dplyr::rename(accession = 'V1', sum.taxonomy = 'V2') %>%
      tidyr::separate(
        .data$sum.taxonomy,
        into = c(
          "superkingdom",
          "phylum",
          "class",
          "order",
          "family",
          "genus",
          "species"
        ),
        sep = ";"
      )
    
    # Convert
    
    ref_db_1_qiime <-
      ref_db_1 %>%  
      dplyr::mutate(
        species = stringr::str_replace_all(.data$species, " ", "_"),
        superkingdom = paste0("k__", .data$superkingdom),
        phylum = paste0("p__", .data$phylum),
        class = paste0("c__", .data$class),
        order = paste0("o__", .data$order),
        family = paste0("f__", .data$family),
        genus = paste0("g__", .data$genus),
        species = paste0("s__", .data$species)
      ) %>%
      tidyr::unite(
        .data$sum.taxonomy,
        c(
          "superkingdom",
          "phylum",
          "class",
          "order",
          "family",
          "genus",
          "species"
        ),
        sep = ";"
      )
    
    taxa_table_path <-
      file.path(out_dir, paste0(ref_db_1_name, "_qiime2_taxonomy.txt"))
    
    utils::write.table(
      ref_db_1_qiime,
      file = taxa_table_path,
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t",
      quote = FALSE
    )
    
  }
