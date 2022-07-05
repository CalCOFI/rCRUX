# A wrapper function for blast_datatable that reads a data table from a path
# Then passes it to blast_datatable
# Then writes the output
# path: a path to a csv from get_blast_seeds

# Value: NULL
rcrux_blast <- function(path, ...) {
    blast_seeds <- read.csv(path)
    blast_datatable(blast_seeds, ...)
    # Write output
    return(NULL)
}

# A function that takes a datatable and returns the results of blasting it
# blast_seeds: a datatable formatted like the output from get_blast_seeds_multi_taxa_or_db

# value: a datatable summarizing the results of blasting each row of the input

# side effects: saves its state in several text files in a particular directory
blast_datatable <- function(blast_seeds, dir, sample_size = 1000, wildcards = "NNNN") {
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
            fasta <- run_blastdbcmd(row)
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
        blastn_output <- run_blastn(aggregate_fasta)
        rbind(output_table, blastn_output)

        # save the state of the blast
        save_state(dir, output_table, unsampled_indices, too_many_Ns, not_in_db)
    }
    return(output_table)
}

# Extracts query parameters from a row and calls blastdbcmd
# row: a row from get_blast_seeds

# Value: a fasta-formatted string
run_blastdbcmd() <- function(row, db_dir) {
    # Extract arguments
    accession <- row$accession
    db <- row$database_used_for_blast
    if(db ==  "refseq_representative_genomes") {
        if(row$superkingdom == "Eukaryota") {
            db <- "/ref_euk_rep_genomes"
        }
        else {
            db <- "/ref_prok_rep_genomes"
        }
    }
    db_path <- paste(db_dir, db, sep = "/")
    db_type <- "nucl"
    forward <- row$forward_stop
    reverse <- row$reverse_stop
    if(forward < reverse) {

    }
    else {
        
    }

    # System call
}

# Runs blastn with the input fasta as a query
# fasta: a fasta-formatted string

# Value: a tibble
run_blastn <- function(fasta) {
    # System call
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

