# Libraries
library(lubridate)
library(XML)
library(httr)
library(ShortRead)
library(primerTree)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ape)
library(tibble)
library(rlist)
library(rlang)
library(taxonomizr)
library(data.table)
library(RCurl)
library(parallel)

##########################
# Main function
# Repeatedly calls run_serial_blast while tracking where it left off

# Calls:
  # make_initial_files
  # load_incomplete_files
  # run_serial_blast

# Apparent global variable dependencies:
  # blast_out
  # to_be_blasted_entries1
  # too_many_Ns
  # too_new_for_you
  # blast_count
  # blasted_number

# Global variables modified/created:
  # None

# Other dependencies:
  # If blast_tracker.txt exists then to_be_blasted_entries.csv must also exist

# Legacy comments:
  # Function to run rCRUX blast on blast_seeds output
  # will also restart where a run left off
  # locked into nucleotide blast and 70% length and ident and e-value....  oops  fix someday but not today

rcrux_blast <- function(file_out_dir, metabarcode_name,
  blast_tools, blast_db, accessionTaxa, number_Ns_in_blast_seed = "NNNN") {
    
  metabarcode_name <- metabarcode_name
  
  # first time running function
  if (!file.exists(paste0(file_out_dir, metabarcode_name, "blast_tracker.txt"))) {
    csv_with_taxonomy_path <- paste0(file_out_dir, metabarcode_name,
                                     "_primerTree_output_with_taxonomy.csv")
    
    # To do: generate an actual error but keep a useful hint
    if (!file.exists(csv_with_taxonomy_path)) {
      print("Run get_blast_seed to run this function")
    }
    
    else {
      print("Starting BLAST")
      make_initial_files(file_out_dir, metabarcode_name, number_Ns_in_blast_seed)
      run_serial_blast(file_out_dir, metabarcode_name, blast_out, to_be_blasted_entries1, too_many_Ns, too_new_for_you, blast_count, blasted_number,blast_tools, blast_db, accessionTaxa, number_Ns_in_blast_seed = "NNNN", dbtype = "'nucl'")
    }
  }
  
  # subsequent times
  else {
    
    # To do: this gets read solely to check its length. Surely there is a better way?
    tbb_tibble <- read_csv(paste0(file_out_dir, metabarcode_name, "_to_be_blasted_entries.csv"))
    
    if (dim(tbb_tibble)[1] == 0){
      print("There is nothing left to BLAST")
    }
    
    else {
      print("resuming BLASTING where you left off")
      load_incomplete_files(file_out_dir, metabarcode_name)
      run_serial_blast(file_out_dir, metabarcode_name, blast_out, to_be_blasted_entries1, too_many_Ns, too_new_for_you, blast_count, blasted_number,blast_tools, blast_db, accessionTaxa, number_Ns_in_blast_seed = "NNNN", dbtype = "'nucl'")
    }  
  }
  
}

##########################
# Creates several files in the directory i_file_out and assigns several global variables
# As the name suggests, it is typically called when those variables and files do not yet exist

# Calls:
  # save_output_as_csv

# Apparent global variable dependencies:
  # None

# Global variables modified/created:
  # blast_out
  # to_be_blasted_entries
  # too_many_Ns
  # too_new_for_you
  # to_be_blasted entries1
    # same as to_be_blasted entries
  # blast_count
    # literally just 1
  # blasted_number
    # literally just 0

# Files written:
  # [prefix]_accessions_with_more_than_[numNs]_in_a_row.csv
  # [prefix]_accessions_not_found_in_your_blast_DB.csv

# Other dependencies:
  # None visible

# Legacy comments:
  # If running blast_rcrux() for the first time need to make the initial files

make_initial_files <- function(i_file_out, i_Metabarcode, number_Ns_in_blast_seed){
  
  colnames <- c("accession",
                "amplicon_length",
                "pident",
                "query_accession",
                "accession_sequence_length", 
                "amplicon_start",
                "amplicon_stop",
                "sequence",
                "evalue",
                "BLAST_db_taxids")
  
  colnames.1 <-c("gi",
                 "accession",
                 "product_length",
                 "mismatch_forward",
                 "mismatch_reverse",
                 "forward_start",
                 "forward_stop",
                 "reverse_start",
                 "reverse_stop",
                 "product_start",
                 "product_stop",
                 "amplicon_length",
                 "taxID",
                 "species",
                 "superkingdom",
                 "kingdom",
                 "phylum",
                 "subphylum",
                 "superclass",
                 "class",
                 "subclass",
                 "order",
                 "family",
                 "subfamily",
                 "genus",
                 "infraorder",
                 "subcohort",
                 "superorder",
                 "superfamily",
                 "tribe",
                 "subspecies",
                 "subgenus",
                 "species group",
                 "parvorder",
                 "varietas")
  
  blast_out <<- data.frame(matrix(ncol = 10, nrow = 0))
  colnames(blast_out) <<- colnames
  
  to_be_blasted_entries <- suppressWarnings(read_csv(paste0(i_file_out, i_Metabarcode, "_primerTree_output_with_taxonomy.csv"))) # figure out error at some point
  # Why as tibble?? It should already by a tibble
  to_be_blasted_entries <<- as_tibble(to_be_blasted_entries)
  
  # Create a global dataframe to record entries with too many Ns
  # Write that dataframe to a csv called [prefix]_accessions_with_more_than_[numNs]_in_a_row.csv
  too_many_Ns <<- data.frame(matrix(ncol = 35, nrow = 0))
  colnames(too_many_Ns) <<- colnames.1
  save_output_as_csv(too_many_Ns, paste0("_accessions_with_more_than_", number_Ns_in_blast_seed ,"_in_a_row"), i_file_out, i_Metabarcode)
  
  # Create a global dataframe to record ???
  # Write that to a csv called [prefix]_accessions_not_found_in_your_blast_DB.csv
  too_new_for_you <<- data.frame(matrix(ncol = 35, nrow = 0))
  colnames(too_new_for_you) <<- colnames.1
  save_output_as_csv(too_new_for_you, "_accessions_not_found_in_your_blast_DB", i_file_out, i_Metabarcode)
  
  to_be_blasted_entries1 <<- to_be_blasted_entries
  
  blast_count <<- 1
  blasted_number <<- 0
  
}

##########################
# Reads several files as tibbles and assigns them to global variables
# Then calls extract_blast_tracker

# Calls:
  # extract_blast_tracker
  
# Apparent global variable dependencies:
  # none

# Global variables modified/created:
  # blast_out
  # to_be_blasted_entries
  # to_be_blasted_entries1
    # same as to_be_blasted_entries
  # too_many_Ns
  # too_new_for_you

# Other dependencies
  # There must be the following files:
    # [prefix]_blasted_entries.csv
    # [prefix]_to_be_blasted_entries.csv
    # [prefix]_accessions_with

# Legacy comments
  # function to load files so blast can pick up where it left off
  # fix duplicated blast out

# number_Ns shouldn't really be a default, this is just a patch to make it work for now until we can do serious refactoring
load_incomplete_files <- function(file_out, Metabarcode, number_Ns_in_blast_seed = "NNNN"){
  blast_out <- suppressWarnings(read_csv(paste0(file_out, Metabarcode, "_blasted_entries.csv")))
  blast_out <<- as_tibble(blast_out)
  
  to_be_blasted_entries <- suppressWarnings(read_csv(paste0(file_out, Metabarcode, "_to_be_blasted_entries.csv"))) 
  to_be_blasted_entries <<- as_tibble(to_be_blasted_entries)
  to_be_blasted_entries1 <<- to_be_blasted_entries
  
  too_many_Ns <- suppressWarnings(read_csv(paste0(file_out, Metabarcode, paste0("_accessions_with_more_than_", number_Ns_in_blast_seed ,"_in_a_row"))))
  too_many_Ns <<- as_tibble(too_many_Ns)
  
  too_new_for_you  <- suppressWarnings(read_csv(paste0(file_out, Metabarcode, "_accessions_not_found_in_your_blast_DB.csv")))
  too_new_for_you <<- as_tibble(too_new_for_you)
  
  extract_blast_tracker(file_out, Metabarcode)
  
}

#################################
# Sample seq_to_blast from to_be_blasted_entries1
# Use blastdbcmd to create a fasta from the sample
# Run blastn with the fasta as an argument
# Saves outputs to too_many_Ns, too_new_for_you, blast_out1, and blast_out
# blast_out1 is results from this blast, blast_out is all of them
# Updates to_be_blasted_entries to reflect things that have been blasted

# To-do
  # Clean up the part that runs blastn
  # See if there is a way to pass a character vector to blastn as the query
  # This can probably be changed to while(nrow(to_be_blasted_entries1) > 0) and remove any logic about when to break

# Calls:
  # save_output_as_csv
  # blast_tracker
  # fetch_rcrux_taxonomy
  # get_fasta_no_hyp
  # run_blastdbcommand
  # run_blastn

# Apparent global variable dependencies:

# Global variables modified/created:
  # None

# Other dependencies:
  # run_blastdbcommand must create a file at the path given by input

# Legacy comments:
  # Function to run serial blast
  # the too many N's and new to you while loop used to be a function.  Fix someday

run_serial_blast <- function(file_out_dir, Metabarcode_name, blast_out, to_be_blasted_entries1, too_many_Ns, too_new_for_you, blast_count, blasted_number,blast_tools, blast_db, accessionTaxa, seq_to_blast = 200, number_Ns_in_blast_seed = "NNNN", dbtype = "nucl"){
  
  blastdbcmd <- paste0(blast_tools,"/bin/blastdbcmd")
  blastn <- paste0(blast_tools,"/bin/blastn")
  
  repeat {
    
    # These variables control whether max_or_seq_to_blast is incremented
    # I'm pretty sure they could go inside the if block below
    i <- 1
    s <- 1

    # Adding these here to make the scope less weird
    max_or_seq_to_blast <- seq_to_blast
    end <- 0
    
    # if there are entries to be blasted
    if (dim(to_be_blasted_entries1)[1] != 0) {
      
      # Check for end condition 
      if (nrow(to_be_blasted_entries1) < seq_to_blast){
        max_or_seq_to_blast <- nrow(to_be_blasted_entries1)
        end <- 1
      }
      
      print(paste0("Get usable accessions for BLAST round number ", blast_count, ":" ))
      
      # Remove the FASTA from the previous round
      unlink(paste0(file_out_dir , "blastdbcmd_test_output.txt"))

      # Variable to track what parts of to_be_blasted_entries1 has been used without modifying it
      unsampled_indices <- c(1:nrow(to_be_blasted_entries1))
      for(i in 1:max_or_seq_to_blast) {
        # Take our sample and extract necessary information
        sample_index <- sample(unsampled_indices, 1)
        unsampled_indices <- unsampled_indices[unsampled_indices!=sample_index]
        sample_row <- to_be_blasted_entries1[sample_index,]
        sample_accession <- sample_row$accession
        sample_db <- paste0(blast_db, "/", sample_row$database_used_for_blast)
        # Correct for the way refseq subdivides data
        if(sample_row$database_used_for_blast ==  "refseq_representative_genomes") {
          if(sample_row$superkingdom == "Eukaryota") {
            sample_db <- paste0(blast_db, "/ref_euk_rep_genomes")
          }
          else {
            sample_db <- paste0(blast_db, "/ref_prok_rep_genomes")
          }
        }

        message(paste0("....trying ", sample_accession))
        blastdbcmd_out_path <- paste0(file_out_dir, "blastdbcmd_test_output_", i, "_.txt")
        run_blastdbcommand(sample_row, blastdbcmd, sample_db, blastdbcmd_out_path, "nucl")

        # Check for problems, record them in the appropriate files
        if (file.info(blastdbcmd_out_path)$size == 0 && dim(to_be_blasted_entries1)[1] != 0) {
          print("....is not in your BLAST database")
          too_new_for_you <- rbind(too_new_for_you, sample_row)
          save_output_as_csv(too_new_for_you, "_accessions_not_found_in_your_blast_DB", file_out_dir, Metabarcode_name)
        }
        else if (any(grepl(number_Ns_in_blast_seed, readLines(blastdbcmd_out_path)))){
          print("....has too Many Ns!")
          too_many_Ns <- rbind(too_many_Ns, subset_for_blast)
          save_output_as_csv(too_many_Ns, paste0("_accessions_with_more_than_", number_Ns_in_blast_seed ,"_in_a_row"), file_out_dir, Metabarcode_name)
          unlink(blastdbcmd_out_path)
        }
      }
    
    # Concatonate blast seed files into a single FASTA
    pattern <- "_.txt$"
    fasta <- readFasta(file_out_dir, pattern)
    writeFasta(fasta, paste0(file_out_dir , "blastdbcmd_test_output.txt"))
    
    # delete files created earlier
    for (j in 1:max_or_seq_to_blast){
      unlink(paste0(file_out_dir , "blastdbcmd_test_output_", j , "_.txt"))
    }
    
    # This variable subtly changes purpose several times through this function
    # The name is admirably concise but not always descriptive
    # I'd like to do some creative refactoring on the uses of this variable
    input <- paste0(file_out_dir, "blastdbcmd_test_output.txt")
    
    # Here's where we do the blastn
    # Could we split this function into two parts? Most things that happen above
      # 1) Are conceptually connected by involving blastdbcmd
      # 2) Seem to write their outputs to one file that can easily be passed to the next function

    # Call blastn with the big FASTA as its query
    if (file.info(input)$size != 0){
      print("....BLASTing")
      print(paste0("....round....", blast_count, ":" ))
      
      blast_out1 <- run_blastn(input, blastn, blast_db)
      unlink(input)
      save_output_as_csv(blast_out1, "_blast_out1", file_out_dir, Metabarcode_name)
      n <- nrow(blast_out1)
      print(paste0("....BLAST ", blast_count, " returned ", n, " rows"))
      
      blast_out <- rbind(blast_out, blast_out1)
      
      # change several tibbles (I think they are tibbles at this point?) to data tables. Why?
      setDT(to_be_blasted_entries1); setDT(blast_out1); setDT(blast_out); setDT(too_many_Ns); setDT(too_new_for_you)
      
      # This appears to be created solely to count the number of rows
      # Perhaps there is a better way to count the number of observations that meet a particular set of criteria?
      # Also we shouldn't be duplicating hits anymore, I think?
      duplicated_blast_out1 <- blast_out1[duplicated(blast_out1$accession) | duplicated(blast_out1$accession, fromLast=TRUE)]
      d <- nrow(duplicated_blast_out1)
      print(paste0("........", d, " duplicated hits"))
      
      # remove duplicates
      blast_out <- setDT(blast_out)[order(-amplicon_length, accession), head(.SD, 1), by = accession]
      save_output_as_csv(blast_out, "_blasted_entries", file_out_dir, Metabarcode_name)

      # not totally sure what this block does, but it appears to log something about the blast
      nt <- nrow(blast_out)
      new <- n - (n + blasted_number  - nt)

      print(paste0("........and ", new, " new blast results"))
      blasted_number <- nt
      print(paste0("...The total number of blast hits for this run is: ", nt))
      
      
      # Remove accessions that match the accessions queries before
      # To-do: review how filtering datatables works in R to better understand this
      to_be_blasted_entries1 <-setDT(to_be_blasted_entries1)[!blast_out, on = .(accession = accession)]
      to_be_blasted_entries1 <-setDT(to_be_blasted_entries1)[!too_many_Ns, on = .(accession = accession)]
      to_be_blasted_entries1 <-setDT(to_be_blasted_entries1)[!too_new_for_you, on = .(accession = accession)]
      
      left <- nrow(to_be_blasted_entries1)
      print(paste0("...There are ", left, " entries that need to be processed"))
      save_output_as_csv(to_be_blasted_entries1, "_to_be_blasted_entries", file_out_dir, Metabarcode_name)
      
      blast_count <- blast_count + 1
      
      
      #update blast_count and blast number to file
      blast_tracker("blast_tracker", file_out_dir, Metabarcode_name, blast_count, blasted_number)
      # Why do you remove it? Shouldn't it be limited in scope?
      rm(subset_for_blast)
      
      if (end == 1) {
        
        fetch_rcrux_taxonomy(blast_out, accessionTaxa, file_out_dir, Metabarcode_name)
        get_fasta_no_hyp(blast_out, file_out_dir, Metabarcode_name)
        
        print(paste0(Metabarcode_name, " BLASTing complete! Find output in: ", file_out_dir, Metabarcode_name))
        break
      }
      
    }
    
    # Perform various cleanup operations
    else {
      fetch_rcrux_taxonomy(blast_out, accessionTaxa, file_out_dir, Metabarcode_name)
      get_fasta_no_hyp(blast_out, file_out_dir, Metabarcode_name)
      print(paste0(Metabarcode_name, " BLASTing complete! Find output in: ", file_out_dir, Metabarcode_name))
      break
    }
  }
  }
}

#################################
# Wraps write.table to point things at the right directory and to append the metabarcode name

save_output_as_csv <- function(file_name, description, file_out, Metabarcode){
  write_to = paste0(file_out, Metabarcode, description, ".csv")
  return(write.table(file_name, file = write_to, row.names=FALSE, sep = ","))
}

#################################
# Updates the values of blast_count and blasted_number based on the contents of blast_tracker.txt

# To-do
  # Move those variables out of the global environment, hold them in the parent function if need be

# Calls:
  # None

# Apparent global variable dependencies:
  # None

# Global variables modified/created:
  # blast_count
  # blasted_number

# Other dependencies:
  # etc_blast_tracker.txt must exist

# Legacy comments (This was grouped with blast_tracker and this comment was above blast_tracker):
  # Function to save current blast run data in .txt files

extract_blast_tracker  <- function(file_out, Metabarcode){
  blast_tracker_file <- file.path(paste0(file_out, Metabarcode, "blast_tracker.txt"))
  get_vars <- file(blast_tracker_file,open="r")
  lines <-readLines(get_vars)
  blast_count <<- as.numeric(lines[1])
  blasted_number <<- as.numeric(lines[2])
  close(get_vars)
}

#################################
# Creates a two line text file that logs two values passed into it
# Used exactly once

# To-do
  # 

# Calls:
  # None

# Apparent global variable dependencies:
  # None

# Global variables modified/created:
  # None

# Other dependencies:
  # 

# Legacy comments:
  # Function to save current blast run data in .txt files

blast_tracker <- function(description, file_out, Metabarcode, blast_count, blasted_number){
  blast_tracker_file <- file.path(paste0(file_out, Metabarcode, description, ".txt"))
  blast_tracker = file(blast_tracker_file,"w");
  write(c(blast_count, blasted_number), blast_tracker ,sep = "/n",append = FALSE, ncolumns = 1);
  close(blast_tracker)
} 

#################################
# Adds taxonomy data to an input datatable and saves the output to a csv

# To-do
  # 

# Calls:
  # get_taxonomizr_from_accession

# Apparent global variable dependencies:
  # None

# Global variables modified/created:
  # None

# Other dependencies:
  # There must be a taxonomizr sql at accessionTaxa

# Legacy comments:
  # Function to add taxonomy to blast output

fetch_rcrux_taxonomy <- function(blast_out, accessionTaxa, file_out_dir, Metabarcode_name){
  Metabarcode_name = Metabarcode_name
  print(paste0("Fetching taxonomy for BLAST output"))
  
  blast_out_with_taxonomy <- get_taxonomizer_from_accession(blast_out, accessionTaxa)
  
  print(paste0("...saving_results"))
  save_output_as_csv(blast_out_with_taxonomy, "_blasted_entries_with_taxonomy", file_out_dir, Metabarcode_name)
}

#################################
# Performs some filtering on an input data table then writes it to a .fasta file

# To-do
  # Review this after reviewing filtering data tables

# Calls:
  # None

# Apparent global variable dependencies:
  # None

# Global variables modified/created:
  # None

# Other dependencies:
  # None

# Legacy comments:
  # None

get_fasta_no_hyp <- function(dupt, file_out_dir, Metabarcode_name){
  dupt_no_hiyp <- dupt %>% mutate(sequence = gsub("-", "", sequence))
  fasta <- character(nrow(dupt_no_hiyp) * 2)
  fasta[c(TRUE, FALSE)] <- paste0(">", dupt_no_hiyp$accession)
  fasta[c(FALSE, TRUE)] <- dupt_no_hiyp$sequence
  writeLines(fasta, paste0(file_out_dir, Metabarcode_name, ".fasta"))
}

#################################
# Given a line from a datatable, the location of blastdbcmd, a blast dbl, and an output path calls blastdbcmd with those arguments
# Uses a helper function to build part of the arguments
# write is for testing purposes

# To-do
  # Right now run_serial_blast does some of the work of parsing information to send here and blastdbcommand_vars does some of the work. Ideally only one function does that work

# Calls:
  # get_blastdbcommand_variables

# Apparent global variable dependencies:
  # None

# Global variables modified/created:
  # None

# Other dependencies:
  # Needs blastdbcmd and a blast database

# Legacy comments:
  # function to run blastdbcommand given the function to get the blastdbcommand variables function above - gets sequence to blast. Saves a file in the working directory called blastdbcmd_test_output.txt:

run_blastdbcommand <- function(data_infile, blastdbcmd, blast_db, output_path, dbtype = "nucl", write = TRUE) {
  if(write) {
    blastdbcommand_vars <- get_blastdbcommand_variables(data_infile)
    # This assumes you have it in your path because that's a cleaner way for this to work for me specifically
    blastdbcmd_out <- system2(command = "blastdbcmd",
                              args = c("-db", blast_db, 
                                      "-dbtype", dbtype,
                                      "-entry",  blastdbcommand_vars,
                                      "-out", output_path),
                              wait = TRUE, stdout = TRUE)
  }
  else {
    blastdbcommand_vars <- get_blastdbcommand_variables(data_infile)
    # This assumes you have it in your path because that's a cleaner way for this to work for me specifically
    blastdbcmd_out <- system2(command = "blastdbcmd",
                              args = c("-db", blast_db, 
                                      "-dbtype", dbtype,
                                      "-entry",  blastdbcommand_vars),
                              wait = TRUE, stdout = TRUE)
  }
  return(blastdbcmd_out)
}

#################################
# Function to run blastn and then parse results as a table
# (I copied the description from the legacy comments)

# To-do
  # Test this complicated blastn command the the CL

# Calls:
  # None

# Apparent global variable dependencies:
  # None

# Global variables modified/created:
  # None

# Other dependencies:
  # The obvious (need to have blastn and have a blastdb)

# Legacy comments:
  # Function to run blastn and then parse results as a table

run_blastn <- function(data_infile , blastn, blast_db, evalue = 1e-6, align = 50000, coverage = 50, perID = 70) {
  blast_db <- paste0(blast_db,"/nt")
  colnames <- c("accession",
                "amplicon_length",
                "pident",
                "query_accession",
                "accession_sequence_length", 
                "amplicon_start",
                "amplicon_stop",
                "sequence",
                "evalue",
                "BLAST_db_taxids") 
  cores <- parallel::detectCores()
  # Assumes in your path because that works better for my specific enviro
  blast_n_output <- system2(command = "blastn", 
                        args = c("-db", blast_db, 
                                 "-query", data_infile, 
                                 "-outfmt", '"6 saccver length pident qacc slen sstart send sseq evalue staxids"', 
                                 "-evalue", evalue,
                                 #"-ungapped", not an inclusive option
                                 "-num_alignments", align,
                                 "-qcov_hsp_perc", coverage,
                                 "-perc_identity", perID,
                                 "-num_threads ", cores),
                        wait = TRUE,
                        stdout = TRUE)
  blast_out1 <- blast_n_output %>%
    as_tibble() %>%
    separate(col = value, into = colnames, sep = "\t", convert = TRUE)
}

#################################
# "Function uses taxonimizer to pull taxonomy from accessions and also collects taxids."
# "Why not just use taxids recovered from blast?  Well blast sometimes pulls multiple taxids."

# To-do
  # Reread after review filtering

# Calls:
  # None

# Apparent global variable dependencies:
  # 

# Global variables modified/created:
  # None

# Other dependencies:
  # taxonomizr sqlite

# Legacy comments:
  # Function uses taxonimizer to pull taxonomy from accessions and also collects taxids.   
  # Why not just use taxids recovered from blast?  Well blast sometimes pulls multiple taxids.  How annoying.... 

get_taxonomizer_from_accession <- function(input, accessionTaxa_path){
  input_taxid <- accessionToTaxa(input$accession, accessionTaxa_path)
  
  input_taxonomy <- getTaxonomy(input_taxid,accessionTaxa_path,desiredTaxa = c("species","superkingdom", "kingdom", "phylum", "subphylum", "superclass", "class", "subclass", "order", "family", "subfamily", "genus", "infraorder", "subcohort", "superorder", "superfamily", "tribe", "subspecies", "subgenus", "species group", "parvorder", "varietas"))
  
  input_taxonomy <- cbind('accession'=input$accession, 'taxID'=input_taxid, input_taxonomy)
  input_taxonomy <- as_tibble(input_taxonomy)
  # Join the blast output and taxonomy tibbles
  return(full_join(input, input_taxonomy, by = "accession"))
}

#################################
# Given a datatable or similar, takes the accession, forward_stop, and reverse_stop
# From the first row and uses them to build a blastdbcmd argument

# To-do
  # There is a weird sharing of responsibilities between this, run_blastdbcommand,
  # and run_serial_blast in processing a row to grab arguments
  # Rework that pipeline

# Calls:
  # None

# Apparent global variable dependencies:
  # 

# Global variables modified/created:
  # 

# Other dependencies:
  # 

# Legacy comments:
  # Function to get the blastdbcommand variables from the datatable.  
  # function chooses accession and range from first row of the infile and formats it for blastdbcommand:

get_blastdbcommand_variables <- function(data_infile) {
  blastdbcommand_accession <- data_infile[1,]$accession
  blastdbcommand_start <- data_infile[1,]$forward_stop
  blastdbcommand_stop <- data_infile[1,]$reverse_stop
  
  if (blastdbcommand_start < blastdbcommand_stop) {
    blastdbcommand_start <- blastdbcommand_start + 1
    blastdbcommand_stop <- blastdbcommand_stop - 1
    return(paste(blastdbcommand_accession," -range ", blastdbcommand_start,"-", blastdbcommand_stop, sep = ''))
  }

  else {
    blastdbcommand_start <- blastdbcommand_start - 1
    blastdbcommand_stop <- blastdbcommand_stop + 1
    return(paste(blastdbcommand_accession," -range ", blastdbcommand_stop,"-", blastdbcommand_start, sep = ''))
  }
}

rcrux_blast("data/home/galoscarleo/blast_seeds_test12S_V5F","12S_V5F",
            "/data/home/galoscarleo/ncbi-blast-2.13.0+/bin",
            "/data/home/galoscarleo/nt",
            "/data/home/galoscarleo/taxonomy/accessionTaxa.sql")