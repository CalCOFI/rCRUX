% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blast_datatable.R
\name{blast_datatable}
\alias{blast_datatable}
\title{Take a datatable and return the results of BLASTing it}
\usage{
blast_datatable(
  blast_seeds,
  save_dir,
  db,
  accession_taxa_path,
  ncbi_bin = NULL,
  force_db = FALSE,
  sample_size = 1000,
  wildcards = "NNNN"
)
}
\arguments{
\item{blast_seeds}{a data.frame formatted like the output from
get_blast_seeds}

\item{save_dir}{a directory in which to create files representing the
current state}

\item{db}{a directory with a blast-formatted database}

\item{accession_taxa_path}{a path to an sql created by
\code{\link[taxonomizr:prepareDatabase]{taxonomizr::prepareDatabase()}}}

\item{ncbi_bin}{the directory that the blast+ suite is in. If NULL, the
program will use your PATH environmental variable to locate them}

\item{force_db}{if true, try to use blast databases that don't appear to
be blast databases}

\item{sample_size}{the number of entries to accumulate into a fasta before
calling blastn}

\item{wildcards}{a character vector representing the number of wildcards to
discard}
}
\value{
A data.frame representing the output of blastn
}
\description{
Given a datatable with the column names of the datatable returned by
\code{\link[RCRUX.dev:get_blast_seeds]{RCRUX.dev::get_blast_seeds()}}, use blastdbcmd to convert entries into
fasta files, then uses blastn to query ncbi databases for those
sequences. It compiles the results of blastn into a data.frame that it
returns. Additionally, it saves its state as text files in a specified
directory with each iteration.
}
\details{
blast_datatable uses blastdbcmd and blastn to find sequences. It samples rows
from \code{blast_seeds} and uses blastdbcmd to find a seqence that corresponds to
the accession number and forward and reverse stops recorded in the table.
blastdbcmd outputs sequences as .fasta-formatted strings, which
blast_datatable concatenates into a multi-line fasta, then passes to blastn
as an argument. blast-datatable repeats this process until no rows remain,
aggregating the results in a single data.frame.
}
\section{Saving data}{
blast_datatable writes intermediate results and metadata about the search to
local files as it goes. This allows the function to resume a partially
completed blast, mitigating the consequences of encountering an
error or experiencing other interruptions. The local files are written to
\code{save_dir} by \code{\link[RCRUX.dev:save_state]{RCRUX.dev::save_state()}}. Manually changing these files is not
suggested as it can change the behavior of blast_datatable. To start from an
incomplete blast_datatable, specify the same save_dir as the incomplete
blast. blast_datable will automatically detect save files and resume from
where it left off.

Warning: If you are resuming from an interrupted blast, make sure you supply
the same data.frame for \code{blast_seeds}. If you intend to start a new blast,
make sure that there is not existing blast save data in the directory you
supply for \code{save_dir}.

Note that blast_datatable does not save intermediate data
from blastdbcmd, so if it is interupted while getting building the fasta to
submit to blastn it will need to repeat some work when resumed. The argument
\code{sample_size} controls the frequency with which it calls blastn, so it can
be used to make blast_datatable save more frequently.
}
