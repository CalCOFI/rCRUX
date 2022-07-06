# A wrapper function for blast_datatable that reads a data table from a path
# Then passes it to blast_datatable
# Then writes the output
# path: a path to a csv from get_blast_seeds
# dir: a directory, passed to blast_datatable for it to save to
# db_dir: a directory containing a blast-formatted database
# expand: logical, determines whether to expand too_many_Ns and not_in_db into real tables

# Value: NULL
rcrux_blast <- function(path, dir, expand = TRUE, ...) {
    blast_seeds <- read.csv(path)
    blast_datatable(blast_seeds, dir, ...)
    # Write output
    return(NULL)
}

# A function that takes a datatable and returns the results of blasting it
# blast_seeds: a datatable formatted like the output from get_blast_seeds_multi_taxa_or_db

# value: a datatable summarizing the results of blasting each row of the input

# side effects: saves its state in several text files in a particular directory
blast_datatable <- function(blast_seeds, dir, db_dir, accession_taxa_path, sample_size = 1000, wildcards = "NNNN") {
    output_table <- data.frame(matrix(ncol = 10, nrow = 0))
    colnames(output_table) <-  c("accession",
                                "amplicon_length",
                                "pident",
                                "query_accession",
                                "accession_sequence_length", 
                                "amplicon_start",
                                "amplicon_stop",
                                "sequence",
                                "evalue",
                                "BLAST_db_taxids")
    num_rounds <- 0
    too_many_Ns <- c()
    not_in_db <- c()
    unsampled_indices <- c(1:nrow(blast_seeds))
    if(file.exists(paste(dir, "unsampled_indices.txt", sep = "/"))) {
        rounds_path <- paste(dir, "num_rounds.txt", sep = "/")
        num_rounds <- as.numeric(readLines(con = rounds_path))

        Ns_path <- paste(dir, "too_many_Ns.txt", sep = "/")
        too_many_Ns <- as.numeric(readLines(con = Ns_path))

        not_in_db_path <- paste(dir, "not_in_db.txt", sep = "/")
        not_in_db <- as.numeric(readLines(con = not_in_db_path))
    }

    while(length(unsampled_indices) > 0) {
        # sample some of them, removing them from the vector
        sample_indices <- smart_sample(unsampled_indices, sample_size)
        unsampled_indices <- unsampled_indices[!(unsampled_indices %in% sample_indices)]

        # run blastdbcmd on each
        # sort results into appropriate buckets
        aggregate_fasta <- ""
        for(index in sample_indices) {
            fasta <- run_blastdbcmd(blast_seeds[index,], db_dir)
            if(nchar(fasta) == 0) {
                append(not_in_db, index)
            }
            else if(grep(wildcards, fasta)) {
                append(too_many_Ns, index)
            }
            else {
                aggregate_fasta <- paste(aggregate_fasta, fasta, sep = "\n")
            }
        }

        # run blastn and aggregate results
        blastn_output <- run_blastn(aggregate_fasta, db_dir)
        output_table <- rbind(output_table, blastn_output)

        # save the state of the blast
        num_rounds <- num_rounds + 1
        save_state(dir, output_table, unsampled_indices, too_many_Ns, not_in_db)
    }
    output_table_taxonomy <- get_taxonomizer_from_accession(output_table, accession_taxa_path)
    return(output_table_taxonomy)
}

# Extracts query parameters from a row and calls blastdbcmd
# row: a row from get_blast_seeds
# ncbi_bin: if not null use it as the parent directory for blastn

# Value: a fasta-formatted string
run_blastdbcmd <- function(query_row, db_dir, ncbi_bin = NULL) {
    # Extract arguments
    accession <- query_row$accession
    db <- query_row$database_used_for_blast
    if(db ==  "refseq_representative_genomes") {
        if(query_row$superkingdom == "Eukaryota") {
            db <- "ref_euk_rep_genomes"
        }
        else {
            db <- "ref_prok_rep_genomes"
        }
    }
    db_path <- paste(db_dir, db, sep = "/")
    db_type <- "nucl"
    forward <- query_row$forward_stop
    reverse <- query_row$reverse_stop
    if(forward < reverse) {
        forward <- forward + 1
        reverse <- reverse -1
    }
    else {
        temp <- forward
        forward <- reverse + 1
        reverse <- temp - 1
    }

    # System call
    if(is.null(ncbi_bin)) {
        fasta <- system2("blastdbcmd", args = c("-db", db_path,
                                                "-dbtype", db_type,
                                                "-entry", accession,
                                                "-range", paste0(forward, "-", reverse)),
                                                stdout = TRUE)
        return(fasta)
    }
    else {
        blastdbcmd <- paste(ncbi_bin, "blastdbcmd", sep = "/")
        fasta <- system2(blastdbcmd, args = c("-db", db_path,
                                                "-dbtype", db_type,
                                                "-entry", accession,
                                                "-range", paste0(forward, "-", reverse)),
                                                stdout = TRUE)
        return(fasta)
    }
}

# Runs blastn with the input fasta as a query
# fasta: a fasta-formatted string
# temp: a file path to write a temporary fasta to
# ncbi_bin: if not null use it as the parent directory for blastn

# Value: a tibble
# Side effects: if temp exists in the working directory, it will be overwritten then deleted
run_blastn <- function(fasta, db_dir, temp = "temp.fasta", ncbi_bin = NULL, evalue = 1e-6, align = 50000, coverage = 50, perID = 70) {
    # This is a hacky workaround to deal with the fact that blastn wants a file path as a query
    # Ideally, we would find a way (perhaps a process substitution?)
    # to pass the string directly to blastn
    # The motivation for writing the fasta here rather than elsewhere is that
    # 1) This way, the script can still generally be visualized functionally
    # 2) It allows for easy changes if we ever figure out an elegant way to do the handoff
    write(fasta, file = temp)

    # Determine arguments
    cores <- parallel::detectCores()

    # System call
    if(is.null(ncbi_bin)) {
        blastn_output <- system2(command = "blastn", 
                                    args = c("-db", blast_db, 
                                    "-query", data_infile, 
                                    "-outfmt", '"6 saccver length pident qacc slen sstart send sseq evalue staxids"', 
                                    "-evalue", evalue,
                                    "-num_alignments", align,
                                    "-qcov_hsp_perc", coverage,
                                    "-perc_identity", perID,
                                    "-num_threads ", cores),
                                    wait = TRUE,
                                    stdout = TRUE)
    }
    else {
        blastn <- paste(ncbi_bin, "blastn", "/")
        blastn_output <- system2(command = blastn, 
                                    args = c("-db", blast_db, 
                                    "-query", data_infile, 
                                    "-outfmt", '"6 saccver length pident qacc slen sstart send sseq evalue staxids"', 
                                    "-evalue", evalue,
                                    "-num_alignments", align,
                                    "-qcov_hsp_perc", coverage,
                                    "-perc_identity", perID,
                                    "-num_threads ", cores),
                                    wait = TRUE,
                                    stdout = TRUE)
    }

    # Format output
    column_names <-  c("accession",
                        "amplicon_length",
                        "pident",
                        "query_accession",
                        "accession_sequence_length", 
                        "amplicon_start",
                        "amplicon_stop",
                        "sequence",
                        "evalue",
                        "BLAST_db_taxids")
    output_table <- blastn_output %>%
                    as_tibble() %>%
                    separate(col = value, into = colnames, sep = "\t", convert = TRUE)
    return(output_table)
}

# Helper function to save the state of a running blast search
# dir: the directory to save in

# Value: NULLL

save_state <- function(dir, output_table, unsampled_indices, too_many_Ns, not_in_db) {
    save_dir <- paste0(dir, "rcrux_blast_data/")
    if(!dir.exists(save_dir)) {
        dir.create(save_dir)
    }
    write.table(output_table,
        file = paste0(save_dir, "output_table.txt"), row.names = FALSE, sep = ",")
    write(unsampled_indices, file = paste0(save_dir, "unsampled_indices.txt"))
    write(too_many_Ns, file = paste0(save_dir, "too_many_Ns.txt"))
    write(not_in_db, file = paste0(save_dir, "not_in_db.txt"))
}

# Wrapper function for sample
# If the sample size is larger than the population, returns the input vector
# Otherwise returns the result of sampling
# x: a vector to sample from
# sized: the size of the sample
# ...: additional arguments to pass to sample

# Value: a vector sampled from x
smart_sample <- function(x, size, ...) {
    if(size > length(x)) {
        return(x)
    }
    else {
        return(sample(x, size, ...))
    }
}

