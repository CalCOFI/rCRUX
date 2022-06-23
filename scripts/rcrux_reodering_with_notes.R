##########################
# Main function
# Repeatedly calls run_serial_blast while tracking where it left off

# Calls:
  # make_initial_files
  # load_incomplete_files
  # run_serial_blast

# Apparent global variable dependencies:
  # blast_out
  # blast_low
  # to_be_blasted_entries
  # too_many_Ns
  # too_new_for_you
  # duplicated_blast_out
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
      run_serial_blast(file_out_dir, metabarcode_name, blast_out, blast_low, to_be_blasted_entries, to_be_blasted_entries1, too_many_Ns, too_new_for_you, duplicated_blast_out, blast_count, blasted_number,blast_tools, blast_db, accessionTaxa, number_Ns_in_blast_seed = "NNNN", dbType = "'nucl'")
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
      run_serial_blast(file_out_dir, metabarcode_name, blast_out, blast_low, to_be_blasted_entries, to_be_blasted_entries1, too_many_Ns, too_new_for_you, duplicated_blast_out, blast_count, blasted_number,blast_tools, blast_db, accessionTaxa, number_Ns_in_blast_seed = "NNNN", dbType = "'nucl'")
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
  # blast_low
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
  
  blast_low <<- data.frame(matrix(ncol = 10, nrow = 0))
  colnames(blast_low) <<- colnames
  
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
#

# Calls:
  # run_blastdbcommand
  # run_blastn
  # save_output_as_csv

# Apparent global variable dependencies:

# Global variables modified/created:

# Other dependencies:
  # run_blastdbcommand must create a file at the path given by input

# Legacy comments:
  # Function to run serial blast
  # the too many N's and new to you while loop used to be a function.  Fix someday

run_serial_blast <- function(file_out_dir, Metabarcode_name, blast_out, blast_low, to_be_blasted_entries, to_be_blasted_entries1, too_many_Ns, too_new_for_you, duplicated_blast_out, blast_count, blasted_number,blast_tools, blast_db, accessionTaxa, seq_to_blast = 200, number_Ns_in_blast_seed = "NNNN", dbType = "'nucl'"){
  
  blastdbcmd <- paste0(blast_tools,"/bin/blastdbcmd")
  blastn <- paste0(blast_tools,"/bin/blastn")
  
  repeat {
    
    # What are these for?
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
      
      # for every number between one and the max number to blast do this loop
      # this is at least one thing s is for
      # Is there a reason for this not to be a for loop?
      while(s <= max_or_seq_to_blast){
        
        # Choose one element of to_be_blasted_entries1 at random
        # Note that this does not alter to_be_blasted_entries1
        subset_for_blast <- sample_n(to_be_blasted_entries1, 1, replace = FALSE)
        
        # The accession number of the sample
        blast_acc <- subset_for_blast$accession[1]
        
        # Determine the database to use
        # This really looks like it could just be a function
        # Default
        blast_db_1 <- paste0(blast_db, "/", subset_for_blast$database_used_for_blast[1])
        # Check for conditions where you should not use the default
        if (subset_for_blast$database_used_for_blast[1] == "refseq_representative_genomes" && subset_for_blast$superkingdom[1] == "Eukaryota") {
          blast_db_1 <- paste0(blast_db, "/ref_euk_rep_genomes")
        }
        else if (subset_for_blast$database_used_for_blast[1] == "refseq_representative_genomes" && subset_for_blast$superkingdom[1] == "Bacteria") {
          blast_db_1 <- paste0(blast_db, "/ref_prok_rep_genomes")
        }
        
        # update the path for blastdbcmd's output, then run it
        print(paste0("....trying ", blast_acc))
        unlink(paste0(file_out_dir , "blastdbcmd_test_output.txt"))
        input <- paste0(file_out_dir, "blastdbcmd_test_output_",s,"_.txt")
        run_blastdbcommand(subset_for_blast, blastdbcmd, blast_db_1, input, dbType)
        
        # for blank blastdbcommand files, check next accession
        # it looks like the truth of the first condition should imply the falsity of the last condition
        # To-do: figure out what on earth this whole block does
        # 
        if ((file.info(input)$size == 0 && dim(to_be_blasted_entries1)[1] != 0) || any(grepl(number_Ns_in_blast_seed, readLines(input)))) {
          # compile observations that are not in the database in too_new_for_you
          if (file.info(input)$size == 0) {
            print("....is not in your BLAST database")
            too_new_for_you <- rbind(too_new_for_you, subset_for_blast)
            # Does this overwrite an existing csv?
            save_output_as_csv(too_new_for_you, "_accessions_not_found_in_your_blast_DB", file_out_dir, Metabarcode_name)
            i = i + 1
          }
          # compile observations with too many Ns in another csv
          # also deletes the blastdbcmd output file then creates an empty one at that location??
          else if (any(grepl(number_Ns_in_blast_seed, readLines(input)))){
            print("....has too Many Ns!")
            too_many_Ns <- rbind(too_many_Ns, subset_for_blast)
            save_output_as_csv(too_many_Ns, paste0("_accessions_with_more_than_", number_Ns_in_blast_seed ,"_in_a_row"), file_out_dir, Metabarcode_name)
            unlink(input)
            file.create(input)
            i = i +1
          }
          # Logically, this should never happen. What is going on??
          else if (file.info(input)$size == 0 && dim(to_be_blasted_entries1)[1] == 0) {
            s = max_or_seq_to_blast
          }
          # if lost a read to too many N's or not in db do not advance the counter, if there are no more hits in the to be blasted file advance the counter one
          # I think the above comment is meant to be associated with the if block below??  
        }

        # To-do: circle back here when you know what i and s do
        if(i == s && s != max_or_seq_to_blast){
          s = s + 1
          i = i + 1
          
        } else if ( i != s && s != max_or_seq_to_blast ) {
          s = s
          i = i - 1
          
        } else {
          s = max_or_seq_to_blast + 1
        }
        
      }  
    }
    
    #concatonate blast seed files
    # Is this meant to be read.FASTA ? If so why is the input formatted wrong?
    # Similarly, is that supposed to be write.FASTA? Perhaps this is because of an outdated package?
    patt <- "_.txt$"
    fasta <- readFasta(file_out_dir, patt)
    writeFasta(fasta, paste0(file_out_dir , "blastdbcmd_test_output.txt"))
    
    # delete files created earlier
    # they either got written to blastdbcmd_test_output.txt
    for (j in 1:max_or_seq_to_blast){
      unlink(paste0(file_out_dir , "blastdbcmd_test_output_", j , "_.txt"))
    }
    
    # This variable subtly changes purpose several times through this function
    # The name is admirably concise but not alwasy descriptive
    # I'd like to do some creative refactoring on the uses of this variable
    input <- paste0(file_out_dir, "blastdbcmd_test_output.txt")
    
    # Here's where we do the blastn
    # Could we split this function into two parts? Most things that happen above
      # 1) Are conceptually connected by involving blastdbcmd
      # 2) Seem to write their outputs to one file that can easily be passed to the next function
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
      duplicated_blast_out1 <- blast_out1[duplicated(blast_out1$accession) | duplicated(blast_out1$accession, fromLast=TRUE)]
      d <- nrow(duplicated_blast_out1)
      print(paste0("........", d, " duplicated hits"))
      
      # remove duplicates
      blast_out <- setDT(blast_out)[order(-amplicon_length, accession), head(.SD, 1), by = accession]
      save_output_as_csv(blast_out, "_blasted_entries", file_out_dir, Metabarcode_name)
      
      # blasted_number appears to be created after being used... potential issue
      # not totally sure what this block does, but it appears to log something about the blast
      nt <- nrow(blast_out)
      new <- n - (n + blasted_number  - nt)

      print(paste0("........and ", new, " new blast results"))
      blasted_number <- nt
      print(paste0("...The total number of blast hits for this run is: ", nt))
      
      
      # remove observations from the primertree output that did not have a blast hit.
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
      
    } else {
      
      fetch_rcrux_taxonomy(blast_out, accessionTaxa, file_out_dir, Metabarcode_name)
      get_fasta_no_hyp(blast_out, file_out_dir, Metabarcode_name)
      print(paste0(Metabarcode_name, " BLASTing complete! Find output in: ", file_out_dir, Metabarcode_name))
      break
    }
  }
  
}