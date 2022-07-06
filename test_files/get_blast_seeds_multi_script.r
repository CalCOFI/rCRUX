# The purpose of this file is to gather the contents of the markdown into an actual R script.

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

# Functions

####################
# Function to download BLAST+.
#
# type works for macox, linux, win64

download_blast_exe <-function(outDir,baseUrl='ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.1/', type='macosx'){
  if (dir.exists(outDir)){
    print(paste0(outDir, " exists"))
  } else {
    dir.create(file.path(outDir))
  }
  
  fileNames<- paste0('ncbi-blast-2.10.1+-x64-',type,'.tar.gz')
  untardir <-  paste0(outDir, 'ncbi-blast-2.10.1+')
  outFiles<-file.path(outDir,fileNames)
  if(dir.exists(untardir)){
    print(paste0(untardir," already exist."))
    return(untardir)
  }
  message('Downloading... This should not take long.')
  urls<-paste(baseUrl,fileNames,sep='/')
  mapply(utils::download.file,urls,outFiles)
  untar(outFiles, exdir=outDir)
  unlink(outFiles)
  return(untardir)
  
}



###########################  
# Function to download BLAST nt database  
# allow for multiple downloads
# make controler to redownload file if there is an error in untar-ing
# make options for multiple blast db downloads.

download_blast_db <- function(blast_db_dir, type='nt'){
  if (dir.exists(blast_db_dir)){
    print(paste0(blast_db_dir, " exists, delete to redownload the BLAST ", type, " directory." ))
  } else {
    dir.create(file.path(blast_db_dir))
    
    require("RCurl") # is this line necessary?
    result <- getURL("ftp://ftp.ncbi.nlm.nih.gov/blast/db/",verbose=TRUE,ftp.use.epsv=FALSE, dirlistonly = TRUE, crlf = TRUE)
    result <- as_tibble((paste("ftp://ftp.ncbi.nlm.nih.gov/blast/db/", strsplit(result, "\r*\n")[[1]], sep = "")))
    result <- rename(result, ftp_name = 1)
    result <- result %>% filter(str_detect(tolower(ftp_name), pattern = paste0("db/", type)))
    result <- result %>% mutate(filenames = ftp_name)
    result <- mutate(result,filenames=sapply(strsplit(result$filenames, split='ftp://ftp.ncbi.nlm.nih.gov/blast/db/', fixed=TRUE),function(x) (x[2])))
    
    message(paste0('Downloading BLAST nt database...This should take a long time.'))
    urls<-paste(result$ftp_name)
    outFiles <- paste0(file.path(blast_db_dir, result$filenames))
    mapply(utils::download.file,urls,outFiles)
    
    files_to_move <- list.files(".", pattern="md5")
    dir.create(file.path(paste0(blast_db_dir,"/md5")))
    for (m in files_to_move){
      file.rename(from = m,  to = paste0(blast_db_dir, "/md5/", m))
    }
    
    files_to_unzip <- list.files(".", pattern="tar.gz")
    for (z in files_to_unzip) {
      print(paste0(blast_db_dir, z))
      untar(paste0(blast_db_dir, z), exdir=blast_db_dir)
      #unlink(z)
    }
  }
  
}





##########################
# Function to save data in .csv files
#

save_output_as_csv <- function(file_name, description, file_out, Metabarcode){
  write_to = paste0(file_out, Metabarcode, description, ".csv")
  return(write.table(file_name, file = write_to, row.names=FALSE, sep = ","))
}



##########################
# Function to save current blast run data in .txt files
#
blast_tracker <- function(description, file_out, Metabarcode, blast_count, blasted_number){
  blast_tracker_file <- file.path(paste0(file_out, Metabarcode, description, ".txt"))
  blast_tracker = file(blast_tracker_file,"w");
  write(c(blast_count, blasted_number), blast_tracker ,sep = "/n",append = FALSE, ncolumns = 1);
  close(blast_tracker)
} 


extract_blast_tracker  <- function(file_out, Metabarcode){
  blast_tracker_file <- file.path(paste0(file_out, Metabarcode, "blast_tracker.txt"))
  get_vars <- file(blast_tracker_file,open="r")
  lines <-readLines(get_vars)
  blast_count <<- as.numeric(lines[1])
  blasted_number <<- as.numeric(lines[2])
  close(get_vars)
}


##########################
# Function to make and save histograms in .csv files
#
Make_hist_save_pdf <- function(infile, description,  file_out_dir, Metabarcode_name){
  pdf(paste0(file_out_dir, Metabarcode_name, description,".pdf"), height = 4, width = 6, onefile=T)
  plot <- hist(infile)
  print(plot)
  dev.off()
}

##########################
# If running blast_rcrux() for the first time need to make the initial files
#

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
  
  #duplicated_blast_out <<- data.frame(matrix(ncol = 10, nrow = 0))
  #colnames(duplicated_blast_out) <<- colnames
  #save_output_as_csv(duplicated_blast_out, "_duplicated_hits_per_blast_run", i_file_out, i_Metabarcode)
  
  
  to_be_blasted_entries <- suppressWarnings(read_csv(paste0(i_file_out, i_Metabarcode, "_primerTree_output_with_taxonomy.csv"))) # figure out error at some point
  to_be_blasted_entries <<- as_tibble(to_be_blasted_entries)
  
  
  too_many_Ns <<- data.frame(matrix(ncol = 35, nrow = 0))
  colnames(too_many_Ns) <<- colnames.1
  save_output_as_csv(too_many_Ns, paste0("_accessions_with_more_than_", number_Ns_in_blast_seed ,"_in_a_row"), i_file_out, i_Metabarcode)
  
  
  too_new_for_you <<- data.frame(matrix(ncol = 35, nrow = 0))
  colnames(too_new_for_you) <<- colnames.1
  save_output_as_csv(too_new_for_you, "_accessions_not_found_in_your_blast_DB", i_file_out, i_Metabarcode)
  
  to_be_blasted_entries1 <<- to_be_blasted_entries
  
  blast_count <<- 1
  blasted_number <<- 0
  #dig <- NULL
  
}

##########################
# function to load files so blast can pick up where it left off
# #fix duplicated blast out

load_incomplete_files <- function(file_out, Metabarcode){
  blast_out <- suppressWarnings(read_csv(paste0(file_out, Metabarcode, "_blasted_entries.csv")))
  blast_out <<- as_tibble(blast_out)
  
  
  #blast_low <- suppressWarnings(read_csv(paste0(file_out, Metabarcode, "_lowest_percent_id_full_length_blast_hits_per_blast_run.csv")))
  #blast_low <<- as_tibble(blast_low)
  
  to_be_blasted_entries <- suppressWarnings(read_csv(paste0(file_out, Metabarcode, "_to_be_blasted_entries.csv"))) 
  to_be_blasted_entries <<- as_tibble(to_be_blasted_entries)
  to_be_blasted_entries1 <<- to_be_blasted_entries
  
  
  too_many_Ns <- suppressWarnings(read_csv(paste0(file_out, Metabarcode, "_accessions_with_more_than_NNNN_in_a_row.csv")))
  too_many_Ns <<- as_tibble(too_many_Ns)
  
  too_new_for_you  <- suppressWarnings(read_csv(paste0(file_out, Metabarcode, "_accessions_not_found_in_your_blast_DB.csv")))
  too_new_for_you <<- as_tibble(too_new_for_you)
  
  #duplicated_blast_out <- suppressWarnings(read_csv(paste0(file_out, Metabarcode, "_duplicated_hits_per_blast_run.csv")))
  #duplicated_blast_out <<- as_tibble(duplicated_blast_out)
  
  extract_blast_tracker(file_out, Metabarcode)
  
}


##########################
# Function to get the blastdbcommand variables from the datatable.  
# function chooses accession and range from first row of the infile and formats it for blastdbcommand:


get_blastdbcommand_variables <- function(data_infile) {
  blastdbcommand_accession <-
    data_infile %>% slice_head( n = 1) %>% pull(accession)
  blastdbcommand_start <-
    data_infile %>% slice_head( n = 1) %>% pull(forward_stop)
  blastdbcommand_stop <-
    data_infile %>% slice_head( n = 1) %>% pull(reverse_stop)
  
  if (blastdbcommand_start < blastdbcommand_stop) {
    blastdbcommand_start <- blastdbcommand_start + 1
    blastdbcommand_stop <- blastdbcommand_stop - 1
    return(paste(blastdbcommand_accession," -range ", blastdbcommand_start,"-", blastdbcommand_stop, sep = ''))
  } else {
    blastdbcommand_start <- blastdbcommand_start - 1
    blastdbcommand_stop <- blastdbcommand_stop + 1
    return(paste(blastdbcommand_accession," -range ", blastdbcommand_stop,"-", blastdbcommand_start, sep = ''))
  }
}


#######################
# function to run blastdbcommand given the function to get the blastdbcommand variables function above - gets sequence to blast. Saves a file in the working directory called blastdbcmd_test_output.txt:

run_blastdbcommand <- function(data_infile, blastdbcmd, blast_db, input, dbType) {
  blastdbcommand_vars <- get_blastdbcommand_variables(data_infile)
  # This assumes you have it in your path because that's a cleaner way for this to work
  blastdbcmd_out <- system2(command = "blastdbcmd",
                            args = c("-db", blast_db, 
                                     "-dbtype", dbType,
                                     "-entry",  blastdbcommand_vars,
                                     "-out", input),
                            wait = TRUE, stdout = TRUE)
}

#######################
#
# function to concatonate blastdbcomd_out files
#

concat_files <- function(file_out_dir){
  wild <- paste0(file_out_dir , "blastdbcmd_test_output_*.txt")
  dataFiles = map(Sys.glob(wild), read_file) 
  sink(paste0(file_out_dir , "blastdbcmd_test_output.txt"))
  writeLines(unlist(lapply(dataFiles, paste, collapse=" ")))
  sink()
  rm(dataFiles)
}




###########################
# Function to check blastdbcommand output. 
# If there are too many Ns -> deletes entry and saves to a specific file
# If output is empty deletes entry
# 

delete_first_line <- function(data_infile, input, too_new_for_you, too_many_Ns, file_out_dir, Metabarcode_name) {
  #print(input)
  if (file.info(input)$size == 0) {
    print("....is not in your BLAST database")
    too_new_for_you1 <<- data_infile %>% slice_head()
    too_new_for_you <<- rbind(too_new_for_you, too_new_for_you1, fill=TRUE) 
    save_output_as_csv(too_new_for_you, "_accessions_not_found_in_your_blast_DB", file_out_dir, Metabarcode_name)
    return(data_infile %>% slice(-1) )# remove top row
  }
  else if (any(grepl(number_Ns_in_blast_seed, readLines(input))) == TRUE) {
    print("....has too Many Ns!")
    too_many_Ns1 <<- data_infile %>% slice_head()
    too_many_Ns <<- rbind(too_many_Ns, too_many_Ns1, fill=TRUE)
    save_output_as_csv(too_many_Ns, paste0("_accessions_with_more_than_", number_Ns_in_blast_seed ,"_in_a_row"), file_out_dir, Metabarcode_name)
    return(data_infile %>% slice(-1) )# remove top row
    unlink(input)
    file.create(input)
  }
  else {
    #print("file has entry")
    return(data_infile)
  }
}

#to_be_blasted_entries1 <- delete_first_line(to_be_blasted_entries, input, too_new_for_you, too_many_Ns, file_out_dir, Metabarcode_name)

############################
# Function to run blastn and then parse results as a table:
#

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
  blast_out1 <- system2(command = "blastn", 
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
                        stdout = TRUE) %>%
    as_tibble() %>% 
    separate(col = value, 
             into = colnames,
             sep = "\t",
             convert = TRUE)
  
}

#arg <- run_blastn( "/Users/limeybean/Documents/GitHub/CRUX_2.0/blastdbcmd_test_output.txt", "/Volumes/20170320_eDNA_seq_data/Blast_and_taxo/Blast_tools/ncbi-blast-2.10.1+/bin/blastn", "/Volumes/20170320_eDNA_seq_data/Blast_and_taxo/Blast_nt_db/nt")

###################
# Function uses taxonimizer to pull taxonomy from accessions and also collects taxids.   
# Why not just use taxids recovered from blast?  Well blast sometimes pulls multiple taxids.  How annoying.... 



get_taxonomizer_from_accession <- function(input, accessionTaxa_path){
  input_taxid <<- accessionToTaxa(input$accession, accessionTaxa_path)
  
  input_taxonomy <<- getTaxonomy(input_taxid,accessionTaxa_path,desiredTaxa = c("species","superkingdom", "kingdom", "phylum", "subphylum", "superclass", "class", "subclass", "order", "family", "subfamily", "genus", "infraorder", "subcohort", "superorder", "superfamily", "tribe", "subspecies", "subgenus", "species group", "parvorder", "varietas"))
  
  input_taxonomy <<- cbind('accession'=input$accession, 'taxID'=input_taxid, input_taxonomy)
  input_taxonomy <<- as_tibble(input_taxonomy)
  # Join the blast output and taxonomy tibbles
  return(full_join(input, input_taxonomy, by = "accession"))
}


# make function for primer_blast



#################################
# Function to run serial blast
# the too many N's and new to you while loop used to be a function.  Fix someday
#

run_serial_blast <- function(file_out_dir, Metabarcode_name, blast_out, blast_low, to_be_blasted_entries, to_be_blasted_entries1, too_many_Ns, too_new_for_you, duplicated_blast_out, blast_count, blasted_number,blast_tools, blast_db, accessionTaxa, seq_to_blast = 200, number_Ns_in_blast_seed = "NNNN", dbType = "nucl"){
  
  blastdbcmd = paste0(blast_tools,"/bin/blastdbcmd")
  blastn = paste0(blast_tools,"/bin/blastn")
  
  repeat {
    
    i = 1
    s = 1
    
    #if there are entries to be blasted... 
    if (dim(to_be_blasted_entries1)[1] != 0) {
      
      #either subset 100 or the number of entries remaining (less than 100)
      if (nrow(to_be_blasted_entries1) >= seq_to_blast){
        max_or_seq_to_blast = seq_to_blast
        end = 0
      } else {
        max_or_seq_to_blast = nrow(to_be_blasted_entries1)
        end = 1
      }
      
      print(paste0("Get usable accessions for BLAST round number ", blast_count, ":" ))
      
      # for every number between one and the max number to blast do this loop
      while(s <= max_or_seq_to_blast){
        
        subset_for_blast <- sample_n(to_be_blasted_entries1, 1, replace = FALSE)
        
        
        
        blast_acc <- subset_for_blast$accession[1]
        
        if (subset_for_blast$database_used_for_blast[1] == "refseq_representative_genomes" && subset_for_blast$superkingdom[1] == "Eukaryota") {
          blast_db_1 <- paste0(blast_db, "/ref_euk_rep_genomes")
        } else if (subset_for_blast$database_used_for_blast[1] == "refseq_representative_genomes" && subset_for_blast$superkingdom[1] == "Bacteria") {
          blast_db_1 <- paste0(blast_db, "/ref_prok_rep_genomes")
        } else {
          blast_db_1 <- paste0(blast_db, "/", subset_for_blast$database_used_for_blast[1])
        }
        
        print(paste0("....trying ", blast_acc))
        unlink(paste0(file_out_dir , "blastdbcmd_test_output.txt"))
        input <- paste0(file_out_dir, "blastdbcmd_test_output_",s,"_.txt")
        suppressWarnings(run_blastdbcommand(subset_for_blast, blastdbcmd, blast_db_1, input, dbType))
        
        # for blank blastdbcommand files, check next accession
        if (file.info(input)$size == 0 && dim(to_be_blasted_entries1)[1] != 0 || any(grepl(number_Ns_in_blast_seed, readLines(input)))){
          #to_be_blasted_entries1 <- delete_first_line(to_be_blasted_entries1, input, too_new_for_you, too_many_Ns, file_out_dir, Metabarcode_name)
          # fix function eventually
          if (file.info(input)$size == 0) {
            print("....is not in your BLAST database")
            too_new_for_you <- rbind(too_new_for_you, subset_for_blast) 
            save_output_as_csv(too_new_for_you, "_accessions_not_found_in_your_blast_DB", file_out_dir, Metabarcode_name)
            i = i +1
          } else if (any(grepl(number_Ns_in_blast_seed, readLines(input))) == TRUE){
            print("....has too Many Ns!")
            too_many_Ns <- rbind(too_many_Ns, subset_for_blast)
            save_output_as_csv(too_many_Ns, paste0("_accessions_with_more_than_", number_Ns_in_blast_seed ,"_in_a_row"), file_out_dir, Metabarcode_name)
            unlink(input)
            file.create(input)
            i = i +1
          } else if (file.info(input)$size == 0 && dim(to_be_blasted_entries1)[1] == 0) {
            
            s = max_or_seq_to_blast
            
          }
          # if lost a read to too many N's or not in db do not advance the counter, if there are no more hits in the to be blasted file advance the counter one
          
          
        }
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
    patt <- "_.txt$"
    # Is this meant to be read.FASTA ? If so why is the input formatted wrong?
    # Similarly, is that supposed to be write.FASTA? Perhaps this is because of an outdated package?
    fasta <- readFasta(file_out_dir, patt)
    writeFasta(fasta, paste0(file_out_dir , "blastdbcmd_test_output.txt"))
    
    
    for (j in 1:max_or_seq_to_blast){
      unlink(paste0(file_out_dir , "blastdbcmd_test_output_", j , "_.txt"))
    }
    
    input <- paste0(file_out_dir, "blastdbcmd_test_output.txt")
    
    if (file.info(input)$size != 0){
      print("....BLASTing")
      print(paste0("....round....", blast_count, ":" ))
      
      blast_out1 <- run_blastn(input, blastn, blast_db)
      unlink(input)
      save_output_as_csv(blast_out1, "_blast_out1", file_out_dir, Metabarcode_name)
      n <- nrow(blast_out1)
      print(paste0("....BLAST ", blast_count, " returned ", n, " rows"))
      
      
      #blast_low1 <- blast_out1 %>% slice_max(amplicon_length) %>% slice_min(pident)
      #blast_low <- rbind(blast_low, blast_low1)
      #save_output_as_csv(blast_low, "_lowest_percent_id_full_length_blast_hits_per_blast_run", file_out_dir, Metabarcode_name)
      
      blast_out <- rbind(blast_out, blast_out1)
      
      
      setDT(to_be_blasted_entries1); setDT(blast_out1); setDT(blast_out); setDT(too_many_Ns); setDT(too_new_for_you)
      
      duplicated_blast_out1 <- blast_out1[duplicated(blast_out1$accession) | duplicated(blast_out1$accession, fromLast=TRUE)]
      d <- nrow(duplicated_blast_out1)
      print(paste0("........", d, " duplicated hits"))
      #duplicated_blast_out <- rbind(duplicated_blast_out, duplicated_blast_out1)
      #save_output_as_csv(duplicated_blast_out, "_duplicated_hits_per_blast_run", file_out_dir, Metabarcode_name)
      
      # remove duplicates
      blast_out <- setDT(blast_out)[order(-amplicon_length, accession), head(.SD, 1), by = accession]
      save_output_as_csv(blast_out, "_blasted_entries", file_out_dir, Metabarcode_name)
      
      nt <- nrow(blast_out)
      new <- n - (n + blasted_number  - nt)
      #return(new)
      print(paste0("........and ", new, " new blast results"))
      blasted_number <- nt
      print(paste0("...The total number of blast hits for this run is: ", nt))
      
      
      # remove observations from the primertree output that did not have a blast hit. 
      to_be_blasted_entries1 <-setDT(to_be_blasted_entries1)[!blast_out, on = .(accession = accession)]
      to_be_blasted_entries1 <-setDT(to_be_blasted_entries1)[!too_many_Ns, on = .(accession = accession)]
      to_be_blasted_entries1 <-setDT(to_be_blasted_entries1)[!too_new_for_you, on = .(accession = accession)]
      
      left <- nrow(to_be_blasted_entries1)
      print(paste0("...There are ", left, " entries that need to be processed"))
      save_output_as_csv(to_be_blasted_entries1, "_to_be_blasted_entries", file_out_dir, Metabarcode_name)
      
      blast_count <- blast_count + 1
      
      
      #update blast_count and blast number to file
      blast_tracker("blast_tracker", file_out_dir, Metabarcode_name, blast_count, blasted_number)
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

get_fasta_no_hyp <- function(dupt, file_out_dir, Metabarcode_name){
  dupt_no_hiyp <- dupt %>% mutate(sequence = gsub("-", "", sequence))
  fasta <- character(nrow(dupt_no_hiyp) * 2)
  fasta[c(TRUE, FALSE)] <- paste0(">", dupt_no_hiyp$accession)
  fasta[c(FALSE, TRUE)] <- dupt_no_hiyp$sequence
  writeLines(fasta, paste0(file_out_dir, Metabarcode_name, ".fasta"))
}

##########################
# Function to run rCRUX blast on blast_seeds output
# will also restart where a run left off
# locked into nucleotide blast and 70% length and ident and e-value....  oops  fix someday but not today

rcrux_blast <- function(file_out_dir, Metabarcode_name, blast_tools, blast_db, accessionTaxa, number_Ns_in_blast_seed = "NNNN"){
  
  Metabarcode_name = Metabarcode_name
  
  if (!file.exists(paste0(file_out_dir, Metabarcode_name, "blast_tracker.txt"))) {
    pto = paste0(file_out_dir, Metabarcode_name, "_primerTree_output_with_taxonomy.csv")
    print(pto)
    if (!file.exists(pto)) {
      print("Run get_blast_seed to run this function")
    } else {
      print("Starting BLAST")
      suppressMessages(make_initial_files(file_out_dir, Metabarcode_name, number_Ns_in_blast_seed))
      print(to_be_blasted_entries1)
      run_serial_blast(file_out_dir, Metabarcode_name, blast_out, blast_low, to_be_blasted_entries, to_be_blasted_entries1, too_many_Ns, too_new_for_you, duplicated_blast_out, blast_count, blasted_number,blast_tools, blast_db, accessionTaxa, number_Ns_in_blast_seed = "NNNN", dbType = "nucl")
      
    }
  } else {
    tbb <- suppressMessages(suppressWarnings(read_csv(paste0(file_out_dir, Metabarcode_name, "_to_be_blasted_entries.csv"))))
    
    #if (file.info(tbb)$size == 0)
    if (dim(tbb)[1] == 0){
      print("There is nothing left to BLAST")
    } else {
      print("resuming BLASTING where you left off")
      suppressMessages(load_incomplete_files(file_out_dir, Metabarcode_name))
      
      run_serial_blast(file_out_dir, Metabarcode_name, blast_out, blast_low, to_be_blasted_entries, to_be_blasted_entries1, too_many_Ns, too_new_for_you, duplicated_blast_out, blast_count, blasted_number,blast_tools, blast_db, accessionTaxa, number_Ns_in_blast_seed = "NNNN", dbType = "nucl")
      
    }  
  }
  
}


######################
# Function to collect BLAST seeds using primerTree
# not sure if it is ready for the big loop wih multiple taxa and primer databases....
#

odds <- function(x) subset(x, x %% 2 != 0)


get_blast_seeds <- function(forward_primer, reverse_primer, file_out_dir, Metabarcode_name, accessionTaxa, num_aligns = 50000, organism = c(organism_to_search), num_permutations = 25, primer_specificity_database = "nt", mismatch = 3, minimum_length = 5, maximum_length = 500){ 
  
  Metabarcode_name = Metabarcode_name
  
  dir.create(file.path(paste0(file_out_dir, Metabarcode_name)))
  
  # search for amplicons using f and r primers
  primer_search_results <- primer_search(forward_primer, reverse_primer, num_aligns = num_aligns, organism =  organism, num_permutations = num_permutations, primer_specificity_database = primer_specificity_database)
  
  # make dataframe
  colnames <- c("gi",
                "accession",
                "product_length",
                "mismatch_forward",
                "mismatch_reverse",
                "forward_start",
                "forward_stop",
                "reverse_start",
                "reverse_stop",
                "product_start",
                "product_stop")
  
  # set up empty tibbles and variables...
  primer_search_blast_out <- data.frame(matrix(ncol = 11, nrow = 0))
  colnames(primer_search_blast_out) <- colnames
  # add break an error messagr -> check primers or use highr taxpnomic rank
  
  url <- list.search(primer_search_results, grepl('http', .), 'character')
  
  t <- length(url)
  times <-odds(1:t)
  
  for (i in times){
    #get url and reformat the output of blast search (there is a smarter way to do this....)
    print(i)
    primer_search_response <- httr::GET(url[[i]])
    print(primer_search_response)
    #parse the blast hits into something human friendly
    primer_search_blast_out1 <- parse_primer_hits(primer_search_response)
    primer_search_blast_out <- rbind(primer_search_blast_out, primer_search_blast_out1) 
    
  }
  
  #make primer_search_blast_out df a tibble
  as_tibble(primer_search_blast_out)
  filter_long_and_short_reads <- primer_search_blast_out %>% filter(mismatch_forward <= mismatch) %>%  filter(mismatch_reverse <= mismatch) %>% filter(product_length >= minimum_length) %>% filter(product_length <= maximum_length) %>% mutate(amplicon_length = product_length - nchar(forward_primer) - nchar(reverse_primer))
  
  # fetch taxonomy associated with the Blast results and arange in alphabetical order starting with species > genus > family > order > class > phylum > superkingdom  - not sure this speeds up blast, but if you are ocd it makes you feel better about life :)
  
  bla <- filter(filter_long_and_short_reads, !grepl(' ', accession))
  
  to_be_blasted_entries <- get_taxonomizer_from_accession(bla, accessionTaxa)
  to_be_blasted_entries <- to_be_blasted_entries %>% arrange(species) %>% arrange(genus) %>% arrange(family)  %>% arrange(order) %>% arrange(class)  %>% arrange(phylum) %>% arrange(superkingdom)
  
  
  
  out <- paste0(file_out_dir, Metabarcode_name, "/")
  
  # save output
  save_output_as_csv(to_be_blasted_entries, "_primerTree_output_with_taxonomy", out, Metabarcode_name)
  Make_hist_save_pdf(primer_search_blast_out$product_length, "_pre_filter_product_lengths_of_primerTree_output",  out, Metabarcode_name)
  save_output_as_csv(primer_search_blast_out, "_raw_primerTree_output", out, Metabarcode_name)
  Make_hist_save_pdf(bla$product_length, "_post_filter_product_lengths_of_primerTree_output",  out, Metabarcode_name)
  
  return(to_be_blasted_entries)
  
  
  
}


#####################################
# Function to add taxonomy to blast output
# 
#

fetch_rcrux_taxonomy <- function(blast_out, accessionTaxa, file_out_dir, Metabarcode_name){
  Metabarcode_name = Metabarcode_name
  print(paste0("Fetching taxonomy for BLAST output"))
  
  blast_out_with_taxonomy <- get_taxonomizer_from_accession(blast_out, accessionTaxa)
  
  print(paste0("...saving_results"))
  save_output_as_csv(blast_out_with_taxonomy, "_blasted_entries_with_taxonomy", file_out_dir, Metabarcode_name)
}



###################
# Function to run a multi taxon multi database search to get primer seeds
#
#

get_blast_seeds_multi_taxa_or_db <- function(forward_primer, reverse_primer, file_out_dir, Metabarcode_name, accessionTaxa, db_to_blast = c("refseq_representative_genomes","nt"),  num_aligns = 50000, organism = c(7777,7742), num_permutations = 25, mismatch = 3, minimum_length = 5, maximum_length = 1000){
  
  #main directory
  dir.create(file.path(paste0(file_out_dir, Metabarcode_name)))
  
  Metabarcode_name = Metabarcode_name
  
  # make main dataframe
  colnames <- c("gi",
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
                "database_used_for_blast",
                "taxid_used_for_blast",
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
  
  # set up empty tibbles and variables...
  all_to_be_blasted_entries <- data.frame(matrix(ncol = 37, nrow = 0))
  colnames(all_to_be_blasted_entries) <- colnames
  
  #loop over db and organism
  dtb = length(db_to_blast)
  org <- length(organism)
  
  for (d in 1:dtb){
    for (o in 1:org){
      print(paste0("primer BLASTing the ", db_to_blast, "database for taxid", organism))
      
      dir.create(file.path(paste0(file_out_dir,Metabarcode_name, "/", Metabarcode_name,"_",db_to_blast[d], "_", organism[o])))
      
      # search for amplicons using f and r primers
      primer_search_results <- primer_search(forward_primer, reverse_primer, num_aligns = num_aligns, organism =  organism[o], num_permutations = num_permutations, primer_specificity_database = db_to_blast[d])
      
      # make dataframe
      colnames <- c("gi",
                    "accession",
                    "product_length",
                    "mismatch_forward",
                    "mismatch_reverse",
                    "forward_start",
                    "forward_stop",
                    "reverse_start",
                    "reverse_stop",
                    "product_start",
                    "product_stop")
      
      # set up empty tibbles and variables...
      primer_search_blast_out <- data.frame(matrix(ncol = 11, nrow = 0))
      colnames(primer_search_blast_out) <- colnames
      # add break an error messagr -> check primers or use highr taxpnomic rank
      
      url <- list.search(primer_search_results, grepl('http', .), 'character')
      
      t <- length(url)
      times <-odds(1:t)
      
      for (i in times){
        #get url and reformat the output of blast search (there is a smarter way to do this....)
        print(i)
        primer_search_response <- httr::GET(url[[i]])
        print(primer_search_response)
        #parse the blast hits into something human friendly
        primer_search_blast_out1 <- parse_primer_hits(primer_search_response)
        primer_search_blast_out <- rbind(primer_search_blast_out, primer_search_blast_out1) 
        
      }
      
      #make primer_search_blast_out df a tibble
      as_tibble(primer_search_blast_out)
      filter_long_and_short_reads <- primer_search_blast_out %>% filter(mismatch_forward <= mismatch) %>%  filter(mismatch_reverse <= mismatch) %>% filter(product_length >= minimum_length) %>% filter(product_length <= maximum_length) %>% mutate(amplicon_length = product_length - nchar(forward_primer) - nchar(reverse_primer)) %>% add_column(database_used_for_blast = db_to_blast[d], taxid_used_for_blast = organism[o])
      
      # fetch taxonomy associated with the Blast results and arange in alphabetical order starting with species > genus > family > order > class > phylum > superkingdom  - not sure this speeds up blast, but if you are ocd it makes you feel better about life :)
      
      bla <- filter(filter_long_and_short_reads, !grepl(' ', accession))
      
      # for testing purposes
      # remove ASAP
      print(bla)
      
      to_be_blasted_entries <- get_taxonomizer_from_accession(bla, accessionTaxa)
      to_be_blasted_entries <- to_be_blasted_entries %>% arrange(species) %>% arrange(genus) %>% arrange(family)  %>% arrange(order) %>% arrange(class)  %>% arrange(phylum) %>% arrange(superkingdom)
      
      
      out <- paste0(file_out_dir, Metabarcode_name, "/", Metabarcode_name,"_",db_to_blast[d], "_", organism[o], "/")
      
      
      # save output
      save_output_as_csv(to_be_blasted_entries, "_primerTree_output_with_taxonomy", out, Metabarcode_name)
      Make_hist_save_pdf(primer_search_blast_out$product_length, "_pre_filter_product_lengths_of_primerTree_output",  out, Metabarcode_name)
      save_output_as_csv(primer_search_blast_out, "_raw_primerTree_output", out, Metabarcode_name)
      Make_hist_save_pdf(bla$product_length, "_post_filter_product_lengths_of_primerTree_output",  out, Metabarcode_name)
      
      #return(to_be_blasted_entries)
      all_to_be_blasted_entries <- rbind(all_to_be_blasted_entries, to_be_blasted_entries)
      
      # delete stuff
      rm(primer_search_results)
      rm(to_be_blasted_entries)
      rm(primer_search_blast_out)
      rm(bla)
      rm(primer_search_blast_out1)
      rm(filter_long_and_short_reads)
      print(paste0("finished and moving on..."))
    }
  }
  
  #print(paste0("Concatonating the primer blast output and removing redundant hits"))
  # remove duplicates
  all_to_be_blasted_entries <- setDT(all_to_be_blasted_entries)[order(-amplicon_length, accession), head(.SD, 1), by = accession]
  out <- paste0(file_out_dir, Metabarcode_name, "/")
  save_output_as_csv(all_to_be_blasted_entries, "_primerTree_output_raw_taxonomy", out, Metabarcode_name)
  all_to_be_blasted_entries_with_taxonomy <- all_to_be_blasted_entries[ taxID != 'NA']
  save_output_as_csv(all_to_be_blasted_entries_with_taxonomy , "_primerTree_output_with_taxonomy", out, Metabarcode_name)
  all_to_be_blasted_entries_too_new_for_you <- all_to_be_blasted_entries[ is.na(taxID)]
  save_output_as_csv(all_to_be_blasted_entries_too_new_for_you, "_primerTree_output_with_taxonomy_not_in_local_db", out, Metabarcode_name)
}


run_serial_blast <- function(file_out_dir, Metabarcode_name, blast_out, blast_low, to_be_blasted_entries, to_be_blasted_entries1, too_many_Ns, too_new_for_you, duplicated_blast_out, blast_count, blasted_number,blast_tools, blast_db, accessionTaxa, seq_to_blast = 1000, number_Ns_in_blast_seed = "NNNN", dbType = "nucl"){
  
  blastdbcmd = paste0(blast_tools,"/bin/blastdbcmd")
  blastn = paste0(blast_tools,"/bin/blastn")
  
  repeat {
    
    i = 1
    s = 1
    
    #if there are entries to be blasted... 
    if (dim(to_be_blasted_entries1)[1] != 0) {
      
      #either subset 100 or the number of entries remaining (less than 100)
      if (nrow(to_be_blasted_entries1) >= seq_to_blast){
        max_or_seq_to_blast = seq_to_blast
        end = 0
      } else {
        max_or_seq_to_blast = nrow(to_be_blasted_entries1)
        end = 1
      }
      
      print(paste0("Get usable accessions for BLAST round number ", blast_count, ":" ))
      
      # for every number between one and the max number to blast do this loop
      while(s <= max_or_seq_to_blast){
        
        subset_for_blast <- sample_n(to_be_blasted_entries1, 1, replace = FALSE)
        
        
        
        blast_acc <- subset_for_blast$accession[1]
        
        if (subset_for_blast$database_used_for_blast[1] == "refseq_representative_genomes" && subset_for_blast$superkingdom[1] == "Eukaryota") {
          blast_db_1 <- paste0(blast_db, "/ref_euk_rep_genomes")
        } else if (subset_for_blast$database_used_for_blast[1] == "refseq_representative_genomes" && subset_for_blast$superkingdom[1] == "Bacteria") {
          blast_db_1 <- paste0(blast_db, "/ref_prok_rep_genomes")
        } else {
          blast_db_1 <- paste0(blast_db, "/", subset_for_blast$database_used_for_blast[1])
        }
        
        
        print(paste0("....trying ", blast_acc))
        unlink(paste0(file_out_dir , "blastdbcmd_test_output.txt"))
        input <- paste0(file_out_dir, "blastdbcmd_test_output_",s,"_.txt")
        suppressWarnings(run_blastdbcommand(subset_for_blast, blastdbcmd, blast_db_1, input, dbType))
        
        # for blank blastdbcommand files, check next accession
        if (file.info(input)$size == 0 && dim(to_be_blasted_entries1)[1] != 0 || any(grepl(number_Ns_in_blast_seed, readLines(input)))){
          #to_be_blasted_entries1 <- delete_first_line(to_be_blasted_entries1, input, too_new_for_you, too_many_Ns, file_out_dir, Metabarcode_name)
          # fix function eventually
          if (file.info(input)$size == 0) {
            print("....is not in your BLAST database")
            too_new_for_you <- rbind(too_new_for_you, subset_for_blast) 
            save_output_as_csv(too_new_for_you, "_accessions_not_found_in_your_blast_DB", file_out_dir, Metabarcode_name)
            i = i +1
          } else if (any(grepl(number_Ns_in_blast_seed, readLines(input))) == TRUE){
            print("....has too Many Ns!")
            too_many_Ns <- rbind(too_many_Ns, subset_for_blast)
            save_output_as_csv(too_many_Ns, paste0("_accessions_with_more_than_", number_Ns_in_blast_seed ,"_in_a_row"), file_out_dir, Metabarcode_name)
            unlink(input)
            file.create(input)
            i = i +1
          } else if (file.info(input)$size == 0 && dim(to_be_blasted_entries1)[1] == 0) {
            
            s = max_or_seq_to_blast
            
          }
          # if lost a read to too many N's or not in db do not advance the counter, if there are no more hits in the to be blasted file advance the counter one
          
          
        }
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
    patt <- "_.txt$"
    fasta <- readFasta(file_out_dir, patt)
    writeFasta(fasta, paste0(file_out_dir , "blastdbcmd_test_output.txt"))
    
    
    for (j in 1:max_or_seq_to_blast){
      unlink(paste0(file_out_dir , "blastdbcmd_test_output_", j , "_.txt"))
    }
    
    input <- paste0(file_out_dir, "blastdbcmd_test_output.txt")
    
    if (file.info(input)$size != 0){
      print("....BLASTing")
      print(paste0("....round....", blast_count, ":" ))
      
      blast_out1 <- run_blastn(input, blastn, blast_db)
      unlink(input)
      save_output_as_csv(blast_out1, "_blast_out1", file_out_dir, Metabarcode_name)
      n <- nrow(blast_out1)
      print(paste0("....BLAST ", blast_count, " returned ", n, " rows"))
      
      
      #blast_low1 <- blast_out1 %>% slice_max(amplicon_length) %>% slice_min(pident)
      #blast_low <- rbind(blast_low, blast_low1)
      #save_output_as_csv(blast_low, "_lowest_percent_id_full_length_blast_hits_per_blast_run", file_out_dir, Metabarcode_name)
      
      blast_out <- rbind(blast_out, blast_out1)
      
      
      setDT(to_be_blasted_entries1); setDT(blast_out1); setDT(blast_out); setDT(too_many_Ns); setDT(too_new_for_you)
      
      duplicated_blast_out1 <- blast_out1[duplicated(blast_out1$accession) | duplicated(blast_out1$accession, fromLast=TRUE)]
      d <- nrow(duplicated_blast_out1)
      print(paste0("........", d, " duplicated hits"))
      #duplicated_blast_out <- rbind(duplicated_blast_out, duplicated_blast_out1)
      #save_output_as_csv(duplicated_blast_out, "_duplicated_hits_per_blast_run", file_out_dir, Metabarcode_name)
      
      # remove duplicates
      blast_out <- setDT(blast_out)[order(-amplicon_length, accession), head(.SD, 1), by = accession]
      save_output_as_csv(blast_out, "_blasted_entries", file_out_dir, Metabarcode_name)
      
      nt <- nrow(blast_out)
      new <- n - (n + blasted_number  - nt)
      #return(new)
      print(paste0("........and ", new, " new blast results"))
      blasted_number <- nt
      print(paste0("...The total number of blast hits for this run is: ", nt))
      
      
      # remove observations from the primertree output that did not have a blast hit. 
      to_be_blasted_entries1 <-setDT(to_be_blasted_entries1)[!blast_out, on = .(accession = accession)]
      to_be_blasted_entries1 <-setDT(to_be_blasted_entries1)[!too_many_Ns, on = .(accession = accession)]
      to_be_blasted_entries1 <-setDT(to_be_blasted_entries1)[!too_new_for_you, on = .(accession = accession)]
      
      left <- nrow(to_be_blasted_entries1)
      print(paste0("...There are ", left, " entries that need to be processed"))
      save_output_as_csv(to_be_blasted_entries1, "_to_be_blasted_entries", file_out_dir, Metabarcode_name)
      
      blast_count <- blast_count + 1
      
      
      #update blast_count and blast number to file
      blast_tracker("blast_tracker", file_out_dir, Metabarcode_name, blast_count, blasted_number)
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

get_fasta_no_hyp <- function(dupt, file_out_dir, Metabarcode_name){
  dupt_no_hiyp <- dupt %>% mutate(sequence = gsub("-", "", sequence))
  fasta <- character(nrow(dupt_no_hiyp) * 2)
  fasta[c(TRUE, FALSE)] <- paste0(">", dupt_no_hiyp$accession)
  fasta[c(FALSE, TRUE)] <- dupt_no_hiyp$sequence
  writeLines(fasta, paste0(file_out_dir, Metabarcode_name, ".fasta"))
}

# Add arguments
get_blast_seeds_multi_taxa_or_db("TAGAACAGGCTCCTCTAG", "TTAGATACCCCACTATGC", "blast_seeds_test", "12S_V5F1",
                                  "/data/home/galoscarleo/taxonomy/accessionTaxa.sql")