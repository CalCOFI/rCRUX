# rCRUX: Generate CRUX metabarcoding reference libraries in R



**Authors:** [Luna Gal](<https://github.com/LunaGal>), [Zachary Gold](<https://github.com/zjgold>), [Ramon Gallego](https://github.com/ramongallego), [Emily Curd](<https://github.com/limey-bean>)<br/>
**Inspiration:**
The late, great [Jesse Gomer](<https://github.com/jessegomer?tab=repositories>). Coding extraordinaire and dear friend.<br/>
**License:**
[GPL-3](https://opensource.org/licenses/GPL-3.0) <br/>
**Support:**
Support for the development of this tool was provided by [CalCOFI](<https://calcofi.org/>), NOAA, Landmark College, and [VBRN](https://vbrn.org/). <br/>
**Acknowledgments:**
This work benefited from the amazing input of many including Lenore Pipes, Sarah Stinson, Gaurav Kandlikar, and Maura Palacios Mejia.


eDNA metabarcoding is increasingly used to survey biological communities using common universal and novel genetic loci. There is a need for an easy to implement computational tool that can generate metabarcoding reference libraries for any locus, and are specific and comprehensive. We have reimagined CRUX ([Curd et al. 2019](https://doi.org/10.1111/2041-210X.13214)) and developed the rCRUX package R system for statistical computing [R Core Team 2021](https://www.r-project.org/) to fit this need by generating taxonomy and fasta files for any user defined locus.  The typical workflow involves using get_seeds_local() or get_seeds_remote() to simulate *in silico* PCR (e.g. [Ye et al. 2012](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-134)) to acquire a set of sequences analogous to PCR products containing metabarcode primer sequences.  The sequences or "seeds" recovered from the *in silico* PCR step are used to search databases for complementary sequence that lack one or both primers. This search step, blast_seeds() is used to iteratively align seed sequences against a local NCBI database for matches using a taxonomic rank based stratified random sampling approach.  This step results in a comprehensive database of primer specific reference barcode sequences from NCBI. Using derep_and_clean_db(), the database is de-replicated by DNA sequence where identical sequences are collapsed into a representative read. If there are multiple possible taxonomic paths for a read, the taxonomic path is collapsed to the lowest taxonomic agreement.


## Typical Workflow
<img src="/flowcharts/rCRUX_overview_flowchart.png" width = 10000 />


## Installation

Install from GitHub:

``` r
# install.packages(devtools)
devtools::install_github("LunaGal/rCRUX")
```

``` r
library(rCRUX)
```

## Dependencies

**NOTE:** These only need to be downloaded once or as NCBI updates databases. rCRUX can access and successfully build metabarcode references using databases stored on external drives. </br>

### BLAST+

NCBI's [BLAST+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/) suite must be locally installed and accessible in the user's path. NCBI provides installation instructions for [Windows](https://www.ncbi.nlm.nih.gov/books/NBK52637/), [Linux](https://www.ncbi.nlm.nih.gov/books/NBK52640/), and [Mac OS](https://www.ncbi.nlm.nih.gov/books/NBK569861/). Version 2.10.1+ is verified compatible with rCRUX.

The following is example script to download blast executables:

```
cd /path/to/Applications

wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.1/ncbi-blast-2.10.1+-x64-macosx.tar.gz

```

This [link](https://community.rstudio.com/t/adding-to-the-path-variable/12066) may help if you are using RStudio and having trouble adding blast+ to your path.


### Blast-formatted database

rCRUX requires a local blast-formatted nucleotide database. These can be user generated or download a pre-formatted database from [NCBI](https://ftp.ncbi.nlm.nih.gov/blast/db/).  NCBI provides a tool (perl script) for downloading databases as part of the blast+ package. A brief help page can be found [here](https://www.ncbi.nlm.nih.gov/books/NBK569850/).

The following shell script can be used to download the blast-formatted nucleotide database.

```

mkdir NCBI_blast_nt

cd NCBI_blast_nt

wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt\*

time for file in *.tar.gz; do tar -zxvf $file; done

cd ..

```

You can test your nt blast database using the following command:

```
blastdbcmd -db '/my/directory/ncbi_nt/nt' -dbtype nucl -entry MN937193.1 -range 499-633
```
If you do not get the following, something went wrong in the build.

```
>MN937193.1:499-633 Jaydia carinatus mitochondrion, complete genome
TTAGATACCCCACTATGCCTAGTCTTAAACCTAGATAGAACCCTACCTATTCTATCCGCCCGGGTACTACGAGCACCAGC
TTAAAACCCAAAGGACTTGGCGGCGCTTCACACCCACCTAGAGGAGCCTGTTCTA
```

***Possible error include but are not limited to:***
1. Partial downloads of database files. Extracting each TAR archive (e.g. nt.00.tar.gz.md5) should result in 8 files with the following extensions(.nhd, .nhi, .nhr, .nin, .nnd, .nni, .nog, and .nsq).  If a few archives fail during download, you can re-download and unpack only those that failed. You do not have to re-download all archives.  

2. You downloaded and built a blast database from ncbi fasta files but did not specify -parse_seqids


The nt database is **~242 GB** (as of 8/31/22) and can take several hours (overnight) to build. Loss of internet connection can lead to partially downloaded files and blastn errors (see above).


### Taxonomizr

rCRUX uses the [taxonomizr](https://cran.r-project.org/web/packages/taxonomizr/vignettes/usage.html) package for taxonomic assignment based on NCBI [Taxonomy id's \(taxids\)](https://www.ncbi.nlm.nih.gov/taxonomy). Many rCRUX functions require a path to a local taxonomizr readable sqlite database. This database can be built using taxonomizr's [prepareDatabase](https://www.rdocumentation.org/packages/taxonomizr/versions/0.8.0/topics/prepareDatabase) function.

This database is **~72 GB** (as of 8/31/22) and can take several hours (overnight) to build. Loss of internet connection can lead to partially downloaded files and taxonomizr run errors.

The following code can be used to build this database:

```
library(taxonomizr)

accession_taxa_sql_path <- "/my/accessionTaxa.sql"
prepareDatabase(accession_taxa_sql_path)

```

**Note:** For poor bandwidth connections, please see the [taxononmizr readme for manual installation](https://cran.r-project.org/web/packages/taxonomizr/readme/README.html) of the accessionTaxa.sql database. If built manually, make sure to delete any files other than the accessionTaxa.sql database (e.g. keeping nucl_gb.accession2taxid.gz leads to a warning message).

# Example pipeline

The following example shows a simple rCRUX pipeline from start to finish. Note that this example will require internet access and considerable database storage (~**314 GB**, see section above), run time (mainly for blastn), and system resources to execute.

**Note:** Blast databases and the taxonomic assignment databases (accessionTaxa.sql) can be stored on external hard drive. It increases run time, but is a good option if computer storage capacity is limited.

There are two options to generate seeds for the database generating blast step blast_seeds_local() or blast_seeds_remote(). The local option is slower, however it is not subject to the memory limitations of using the NCBI primer_blast API. The local option is recommended if the user is building a large database, wants to include any [taxid](https://www.ncbi.nlm.nih.gov/taxonomy) in the search, and has many degenerate sites in their primer set.

## get_seeds_local()

This example uses default parameters, with the exception of evalue to minimize run time.

```

forward_primer_seq = "TAGAACAGGCTCCTCTAG"

reverse_primer_seq =  "TTAGATACCCCACTATGC"

output_directory_path <- "/my/directory/12S_V5F1_local_111122_e300" # path to desired output directory

metabarcode_name <- "12S_V5F1" # desired name of metabarcode locus

accession_taxa_sql_path <- "/my/directory/accessionTaxa.sql" # path to taxonomizr sql database

blast_db_path <- "/my/directory/ncbi_nt/nt"  # path to blast formatted database


get_seeds_local(forward_primer_seq,
                 reverse_primer_seq,
                 output_directory_path,
                 metabarcode_name,
                 accession_taxa_sql_path,
                 blast_db_path, evalue = 300)

```

Two output .csv files are automatically created at this path based on the arguments passed to get_seeds_local.  One includes all unfiltered output the other is filtered based on user defined parameters and includes taxonomy.

A unique taxonomic rank summary file is also generated (e.g. the number of unique phyla, class, etc in the blast hits). If a taxonomic rank category contains NA's, they will be counted as a single unique rank. Sequence availability in NCBI for a given taxid is a limiting factor.

Example output can be found [here](/examples/12S_V5F1_generated_11-11-22).


**If BLAST+ is not in your path do the following:**


```

get_seeds_local(forward_primer_seq,
                 reverse_primer_seq,
                 output_directory_path,
                 metabarcode_name,
                 accession_taxa_sql_path,
                 blast_db_path, evalue = 300,
                 ncbi_bin = "/my/directory/blast+_folder/bin")

```


## get_seeds_remote()

This example uses default parameters to minimize run time.

Searching jawless vertebrates (taxid: "1476529") and jawed vertebrates (taxid: "7776").

```

forward_primer_seq = "TAGAACAGGCTCCTCTAG"

reverse_primer_seq =  "TTAGATACCCCACTATGC"

output_directory_path <- "/my/directory/12S_V5F1_remote_111122" # path to desired output directory

metabarcode_name <- "12S_V5F1" # desired name of metabarcode locus

accession_taxa_sql_path <- "/my/directory/accessionTaxa.sql" # path to taxonomizr sql database



get_seeds_remote(forward_primer_seq,
          reverse_primer_seq,
          output_directory_path,
          metabarcode_name,
          accession_taxa_sql_path,
          organism = c("1476529", "7776"),
          return_table = FALSE)

```

**Note:** When using default parameters only 1047 hits are returned from NCBI's primer blast (run 11-11-22). Returns hit sizes and contents are variable depending on parameters, random blast sampling, and database updates.

Two output .csv files are automatically created at this path based on the arguments passed to get_seeds_remote.  One includes all unfiltered output the other is filtered based on user defined parameters and includes taxonomy.

A unique taxonomic rank summary file is also generated (e.g. the number of unique superkingdon, phyla, class, etc in the blast hits). If a taxonomic rank category contains NA's, they will be counted as a single unique rank.

Sequence availability in NCBI for a given taxid is a limiting factor, as are degenerate bases and API memory allocation.

[Modifying defaults can increase the number of returns by orders of magnitude.](#get_seeds_local)

Example output can be found [here](/examples/12S_V5F1_generated_11-11-22).

## blast_seeds()


Iterative searches are based on a stratified random sampling unique taxonomic groups for a given rank from the get_seeds_local or get_seeds_remote output table. For example, the default is to randomly sample one read from each genus.  The user can select any taxonomic rank present in the get_seeds_local output table. The number of seeds selected may cause blastn to exceed the users available RAM, and for that reason the user can choose the maximum number of reads to blast at one time (max_to_blast, default = 1000). blast_seeds will subsample each set of seeds based on max_to_blast and process all seeds before starting a new search for seeds to blast. It saves the output from each round of blastn.  


```

seeds_output_path <- '/my/directory/12S_V5F1_remote_111122/12S_V5F1_filtered_get_seeds_remote_output_with_taxonomy.csv' # this is output from get_seeds_local or get_seeds_remote

blast_db_path <- "/my/directory/blast_database/nt"  # path to blast formatted database

accession_taxa_sql_path <- "/my/directory/accessionTaxa.sql"  # path to taxonomizr sql database

output_directory_path <- '/my/directory/12S_V5F1_remote_111122/' # path to desired output directory

metabarcode_name <- "12S_V5F1"  # desired name of metabarcode locus


blast_seeds(seeds_output_path,
            blast_db_path,
            accession_taxa_sql_path,
            output_directory_path,
            metabarcode_name)    

```


**Note:** After each round of blast, the system state is saved. If the script is terminated after a full round of blast, the user can pick up where they left off. The user can also change parameters at this point (e.g. change the max_to_blast or rank)

The output includes a summary table of all unique blast hits (summary.csv), a multi fasta file of all unique hits (metabarcode_name_.fasta), a taxonomy file of all unique hits (metabarcode_name_taxonomy.txt), a unique taxonomic rank summary file (metabarcode_name_taxonomic_rank_counts.txt), a list of all of the accessions not present in your blast database (e.g. relevant if you ran get_seeds_remote; blastdbcmd_failed.csv), and a list of accessions with 4 or more Ns in a row (default for that parameter is wildcards = "NNNN"; too_many_ns.csv). The default number of reads to blast per rank is 1 (default for that parameter is sample_size = 1). The script will error out if the user asks for more reads per rank than exist in the blast seeds table.     



**If BLAST+ is not in your path do the following**


```

blast_seeds(seeds_output_path,
            blast_db_path,
            accession_taxa_sql_path,
            output_directory_path,
            metabarcode_name,
            ncbi_bin = "/my/directory/blast+_folder/bin")

```

Example output can be found [here](/examples/12S_V5F1_generated_11-11-22).

**Note:** There will be variability between runs due to primer blast return parameters and random sampling of the blast seeds table that occurs during blast_seeds. However, variability can be decreased by changing parameters (e.g. randomly sampling species rather than genus will decrease run to run variability).

## derep_and_clean_db()


This function takes the output of blast_seeds and de-replicates identical sequences and collapses ambiguous taxonomy to generate a clean reference database.

```

output_directory_path <- '/my/directory/12S_V5F1_remote_111122/' # path to desired output directory

summary_path <- "/my/directory/12S_V5F1_remote_111122/blast_seeds_output/summary.csv" # this is the path to the output from blast_seeds

derep_and_clean_db(output_directory_path, summary_path, metabarcode_name)


```

**Note:** Accessions with the same sequence are collapsed into a representative sequence. If those accessions have different taxids (taxonomic paths), we determine the lowest taxonomic agreement across the multiple accessions with an identical sequence. For example, for the MiFish 12S locus, nearly all rockfishes in the genus *Sebastes* have identical sequences. Instead of including ~110 identical reference sequences, one for each individual species, we report a single representative sequence with a lowest common taxonomic agreement of the genus *Sebastes*. This prevents classification bias for taxa with more sequences and also provides accurate taxonomic resolution within the reference database.

We exclude all sequences with taxids that are NA. Such sequences are not immediately useful for classification of metabarcoding sequences. However, we caution that such results can be indicative of off target amplification of a given primer set. For example, the MiFish 12S primer set amplifies uncultured marine bacteria among other taxa (taxid = NA) indicating off target amplification of non-fish taxa. These sequences are saved in the References_with_NA_for_taxonomic_ranks.csv file.

The result of this function is a final clean reference database file set composed of a paired metabarcode_name_derep_and_clean.fasta and metabarcode_name_derep_and_clean_taxonomy.txt. A summary file of the number of unique taxonomic ranks is also generated: metabarcode_name_derep_and_clean_unique_taxonomic_rank_counts.txt. In addition, all representative sequences and associated accessions are saved in Sequences_with_lowest_common_taxonomic_path_agreement.csv, Sequences_with_mostly_NA_taxonomic_paths.csv,
Sequences_with_multiple_taxonomic_paths.csv, and
Sequences_with_single_taxonomic_path.csv files. These files allow for the traceback of representative sequences to multiple accessions.


Example output can be found [here](/examples/12S_V5F1_generated_11-11-22).


# Detailed Explanation For The Major Functions

## [get_seeds_local](https://limey-bean.github.io/get_seeds_local)

<img src="/flowcharts/get_seeds_local-flowchart.png" width = 10000 />

### Overview
[get_seeds_local](https://limey-bean.github.io/get_seeds_local) takes a set of forward and reverse primer sequences and generates .csv summaries of data returned from a locally run adaptation of [NCBI's primer blast](https://www.ncbi.nlm.nih.gov/tools/primer-blast/). This function performs like *in silicon* to find possible full length barcode sequences containing forward and reverse primer matches. It also generates a count of unique instances of taxonomic ranks (Phylum, Class, Order, Family, Genus, and Species) found in the output.

This script is a local interpretation of [get_seeds_remote](https://limey-bean.github.io/get_seeds_remote) that avoids querying NCBI's [primer BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/) tool. Although it is slower than remotely generating blast seeds, it is not subject to the arbitrary throttling of jobs that require significant memory.

### Expected Output
It creates a `get_seeds_local` directory at `output_directory_path` if one doesn't yet exist, then creates a subdirectory inside `output_directory_path` named after `metabarcode_name`. It creates three files inside that directory. One represents the unfiltered output and another represents the output after filtering with user modifiable parameters and with appended taxonomy. Also generated is a summary of unique taxonomic ranks after filtering.

### Detailed Steps
get_seeds_local passes the forward and reverse primer sequence for a given PCR product to [run_primer_blastn](https://limey-bean.github.io/run_primer_blastn). run_primer_blastn takes a fasta file containing only the forward and reverse primer, and individually queries them against a blast formatted database using the task "blastn_short". The result is an output table with the following columns of data: qseqid (query subject id), sgi (subject gi), saccver (subject accession version), mismatch (number of mismatches between the subject a query), sstart (subject start), send (subject end), staxids (subject taxids).

The returned blast hits for each sequence are matched and checked to see if they generate plausible amplicon (e.g. amplify the same accession and are in the correct orientation to produce a PCR product). These hits are written to a file with the suffix `_unfiltered_get_seeds_local_output.csv`.  These hits are further filtered for length and number of mismatches.

Taxonomy is appended to these filtered hits using [get_taxonomizr_from_accession](https://limey-bean.github.io/get_taxonomizr_from_accession). The results are written to to file with the suffix `_filtered_get_seeds_local_output_with_taxonomy.csv`. The number of unique instances for each rank in the taxonomic path for the filtered hits are tallied (NAs are counted once per rank) and written to a file with the suffix `_filtered_get_seeds_local_unique_taxonomic_rank_counts.txt`.

**Note:**
Information about the blastn parameters can be found in run_primer_blast, and by accessing blastn -help in your terminal.  Default parameters were optimized to provide results similar to those generated through remote blast via primer-blast as implemented in [iterative_primer_search](https://limey-bean.github.io/iterative_primer_search) and [modifiedPrimerTree_Functions](https://limey-bean.github.io/modifiedPrimerTree_Functions). The number of alignments returned for a given blast search is hardcoded at "-num_alignments", "10000000",


### Parameters
**forward_primer_seq**
+ passed to primer_to_fasta, which turns it into fasta
        file to be past to get_seeds_local.
+       e.g. forward_primer_seq <- "TAGAACAGGCTCCTCTAG"
**reverse_primer_seq**
+ passed to primer_to_fasta, which turns it into fasta
        file to be past to get_seeds_local.
+       e.g. reverse_primer_seq <-  "TTAGATACCCCACTATGC"
**output_directory_path**
+ the parent directory to place the data in.
+       e.g. "/path/to/output/12S_V5F1_local_111122_e300_111122"
**metabarcode_name**
+ used to name the subdirectory and the files. get_seeds_local appends
        metabarcode_name to the beginning of each of the files it
        generates.
+       e.g. metabarcode_name <- "12S_V5F1"
**accession_taxa_sql_path**
+ the path to sql database created by taxonomizr
+       e.g. accession_taxa_sql_path <- "/my/accessionTaxa.sql"
**mismatch**
+ the highest acceptable mismatch value per hit. get_seeds_local removes each
        row with a mismatch greater than the specified value.
+       The default is mismatch = 6
**minimum_length**
+ get_seeds_local removes each row that has a value less than
        minimum_length in the product_length column.
+       The default is minimum_length = 5
**maximum_length**
+ get_seeds_local removes each row that has a
        value greater than maximum_length in the product_length column.
+       The default is maximum_length = 500
**blast_db_path**
+ the path to a directory with a blast-formatted database.
+       e.g blast_db_path <- "/my/ncbi_nt/nt"
**task**
+ passed to [run_primer_blastn](https://limey-bean.github.io/run_primer_blastn) the task for blastn to
        perform
+       The default is "blastn_short" - which is optimized for searches with queries < 50 bp
**word_size**
+ passed to [run_primer_blastn](https://limey-bean.github.io/run_primer_blastn) is the fragment size
        used for blastn search - smaller word sizes increase sensitivity and
        time of the search.
+       The default is word_size =  7
**evalue**
+ passed to [run_primer_blastn](https://limey-bean.github.io/run_primer_blastn) is the number of
        expected hits with a similar quality score found by chance.
+       The default is evalue = 3e-7
**coverage**
+ passed to [run_primer_blastn](https://limey-bean.github.io/run_primer_blastn) is the minimum
        percent of the query length recovered in the subject hits.
+       The default is coverage = 90
**perID**
+ passed to [run_primer_blastn](https://limey-bean.github.io/run_primer_blastn) is the minimum percent
        identity of the query relative to the subject hits.
+       The default is perID = 50
**reward**
+ passed to [run_primer_blastn](https://limey-bean.github.io/run_primer_blastn) is the reward for
        nucleotide match.
+       The default is reward = 2
**ncbi_bin**
+ passed to [run_primer_blastn](https://limey-bean.github.io/run_primer_blastn) is the path to blast+
        tools if not in the user's path. Specify only if blastn and is not in
        your path.
+       The default is ncbi_bin = NULL - if not specified in path do the following: ncbi_bin = "/my/local/blast+_folder/bin".

### Example

```

forward_primer_seq = "TAGAACAGGCTCCTCTAG"
reverse_primer_seq =  "TTAGATACCCCACTATGC"
output_directory_path <- "/my/directory/12S_V5F1_local_111122_species_750"
metabarcode_name <- "12S_V5F1" # desired name of metabarcode locus
accession_taxa_sql_path <- "/my/directory/accessionTaxa.sql"
blast_db_path <- "/my/directory/ncbi_nt/nt"

get_seeds_local(forward_primer_seq,
                reverse_primer_seq,
                output_directory_path,
                metabarcode_name,
                accession_taxa_sql_path,
                blast_db_path,
                minimum_length = 80,
                maximum_length = 150)

# adjusting the minimum_length and maximum_length parameters reduces the number of total hits by removing reads that could result from off target amplification

```


## [get_seeds_remote](https://limey-bean.github.io/get_seeds_remote)

<img src="/flowcharts/get_seed_remote-flowchart.png" width = 10000 />

### Overview

[get_seeds_remote](https://limey-bean.github.io/get_seeds_remote) takes a set of forward and reverse primer sequences and generates .csv summaries of [NCBI's primer blast](https://www.ncbi.nlm.nih.gov/tools/primer-blast/) data returns. Only full length barcode sequences containing primer matches are captured. It also generates a count of unique instances of taxonomic ranks (Phylum, Class, Order, Family, Genus, and Species) captured in the seed library.

This script uses [iterative_primer_search](https://limey-bean.github.io/iterative_primer_search) to perform tasks. Its parameters are very similar to primerTree's primer_search(), but it takes vectors for organism and for database and performs a primer search for each combination. For each combination it calls [modifiedPrimerTree_Functions](https://limey-bean.github.io/modifiedPrimerTree_Functions), which is a modified versions of primerTree's primer_search() and primerTree's parse_primer, to query NCBI's [primer BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/) tool, filters the results, and aggregates them into a single data.frame.

It downgrades errors from primer_search and parse_primer_hits into warnings. This is useful when searching for a large number of different combinations, allowing the function to output successful results.

### Expected Output
It creates a directory `get_seeds_remote` in the `output_directory_path`. It creates three files inside that directory. One represents the unfiltered output and another represents the output after filtering with user modifiable parameters and with appended taxonomy. Also generated is a summary of unique taxonomic ranks after filtering.


### Detailed Steps

get_seeds_remote passes the forward and reverse primer sequence for a given
PCR product to [iterative_primer_search](https://limey-bean.github.io/iterative_primer_search) along with the taxid(s) of
the organism(s) to blast, the database to search, and many additional possible
parameters to NCBI's primer blast tool (see Note below). Degenerate primers
are converted into all possible non degenerate sets and a user defined maximum
number of primer combinations is passed to to the API using [modifiedPrimerTree_Functions](https://limey-bean.github.io/modifiedPrimerTree_Functions). Multiple taxids are searched independently, as are multiple databases (e.g. c('nt', 'refseq_representative_genomes'). The data are parsed and stored in a dataframe, which is also written to a file with the suffix `_unfiltered_get_seeds_remote_output.csv`.

These hits are further filtered using [filter_primer_hits](https://limey-bean.github.io/filter_primer_hits) to
calculate and append amplicon size to the dataframe. Only hits that pass with default
or user modified length and number of mismatches parameters are retained.

Taxonomy is appended to these filtered hits using
[get_taxonomizr_from_accession](https://limey-bean.github.io/get_taxonomizr_from_accession). The results are written to
to file with the suffix `_filtered_get_seeds_remote_output_with_taxonomy.csv`.
The number of unique instances for each rank in the taxonomic path for the
filtered hits are tallied (NAs are counted once per rank) and written to a
file with the suffix `_filtered_get_seeds_local_remote_taxonomic_rank_counts.txt`

**Notes:**
get_seeds_remote passes many parameters to NCBI's primer blast tool. See below for more information.

primer BLAST defaults to homo sapiens, so it is important that you supply a specific organism or organisms. NCBI's taxids can be found [here](https://www.ncbi.nlm.nih.gov/taxonomy). You can specify multiple organism by passing a character vector containing each of the options, like in the example below.

Often NCBI API will throttle higher taxonomic ranks (Domain, Phylum, etc.). One work around is to supply multiple lower level taxonomic ranks (Class, Family level, etc.) or use get_seeds_local.

### Parameters

**forward_primer_seq**
+ passed to primer_search, which turns it into a list of
        each primer it could be based on its degenerate primers, then passes
        each one in turn to NCBI.
+       e.g. forward_primer_seq <- "TAGAACAGGCTCCTCTAG"
**reverse_primer_seq**
+ passed to primer_search, which turns it into a list of
        each primer it could be based on its degenerate primers, then passes
        each one in turn to NCBI.
+        e.g. reverse_primer_seq <-  "TTAGATACCCCACTATGC"
**output_directory_path**
+ the parent directory to place the data in.
+        e.g. "/path/to/output/12S_V5F1_remote_111122"
**metabarcode_name**
+ used to name output files. get_seeds_remote appends
        metabarcode_name to the beginning of each of the two files it generates.
+       e.g. metabarcode_name <- "12S_V5F1"
**accession_taxa_sql_path**
+ the path to sql created by taxonomizr
+         e.g. accession_taxa_sql_path <- "/my/accessionTaxa.sql"
**organism**
+ a vector of character vectors. Each character vector is
        passed in turn to primer_search, which passes them to NCBI.
        get_seeds_remote aggregates all of the results into a single file
+       e.g. organism = c("1476529", "7776"))
+       Note: increasing taxonomic rank (e.g. increasing from order to class) for this parameter can maximize primer hits, but can also lead to API run throttling due to memory limitations
**num_permutations**
+ the number of primer permutations to search, if the
        degenerate bases cause more than this number of permutations to exist,
        this number will be sampled from all possible permutations.
+       The default is num_permutations = 50  
+       Note for very degenerate bases, searches may be empty due to poor mutual matches for a given forward and reverse primer combination.
**mismatch**
+ the highest acceptable mismatch value. parse_primer_hits
        returns a table with a mismatch column. get_seeds_remote removes each
        row with a mismatch greater than the specified value.
+       The default is mismatch = 3
+       Note: this is smaller than get_seeds_local because of differences in mismatch calculation between function.
**minimum_length**
+ parse_primer_hits returns a table with a product_length
        column. get_seeds_remote removes each row that has a value less than
        minimum_length in the product_length column.
+       The default is minimum_length = 5
**maximum_length**
+ parse_primer_hits returns a table with a
        product_length column. get_seeds_remote removes each row that has a
        value greater than maximum_length in the product_length column
+       The default is maximum_length = 500
**primer_specificity_database**
+ passed to primer_search, which passes it to NCBI.  
+       The default is primer_specificity_database = 'nt'.
**HITSIZE**
+ a primer BLAST search parameter. Set to a high vlaue to maximize the
        number of observations returned.
+       The default HITSIZE = 50000
+       Note: increasing this parameter can maximize primer hits, but can also lead to API run throttling due to memory limitations
**NUM_TARGETS_WITH_PRIMERS**
+ a primer BLAST search parameter set high to
        maximize the number of observations returned.
+       The default is NCBI NUM_TARGETS_WITH_PRIMERS = 1000
+       Note: increasing this parameter can maximize primer hits, but can also lead to API run throttling due to memory limitations
**...**
+ additional arguments passed to [modifiedPrimerTree_Functions](https://limey-bean.github.io/modifiedPrimerTree_Functions). See [NCBI primer-blast tool]](https://www.ncbi.nlm.nih.gov/tools/primer-blast/) for more information.


**Check NCBI's primer blast for additional search options**

get_seeds_remote passes many parameters to NCBI's primer blast tool. You can match the parameters to the fields available in the GUI here. First, use your browser to view the page source. Search for the field you are interested in by searching for the title of the field. It should be enclosed in a tag. Inside the label tag, it says for = "<name_of_parameter>". Copy the string after for = and add it to get_seeds_remote as the name of a parameter, setting it equal to whatever you like.

As of 2022-08-16, the primer blast GUI contains some options that are not implemented by primer_search. The [table below] documents some of the available options.

| Name                                   |       Default  |
|----------------------------------------|----------------|
| PRIMER_SPECIFICITY_DATABASE            | nt             |
| EXCLUDE_ENV                            | unchecked      |
| ORGANISM                               | Homo sapiens   |
| TOTAL_PRIMER_SPECIFICITY_MISMATCH      | 1              |
| PRIMER_3END_SPECIFICITY_MISMATCH       | 1              |
| TOTAL_MISMATCH_IGNORE                  | 6              |
| MAX_TARGET_SIZE                        | 4000           |
| HITSIZE                                | 50000          |
| EVALUE                                 | 30000          |
| WORD_SIZE                              | 7              |
| NUM_TARGETS_WITH_PRIMERS               | 1000           |
| MAX_TARGET_PER_TEMPLATE                | 100            |



You can check [primerblast](https://www.ncbi.nlm.nih.gov/tools/primer-blast/) for more information on how to modify search options. For example, if want you to generate a larger hitsize, open the source of the primer designing tool and look for that string. You find the following:

```
<label for="HITSIZE" class="m ">Max number of sequences returned by Blast</label>
         <div class="input ">
                      <span class="sel si">
                      <select name="HITSIZE" id="HITSIZE" class= "opts checkDef" defVal="50000" >
                        <option  value="10">10</option>
                        <option  value="50">50</option>
                        <option  value="100">100</option>
                        <option  value ="250">250</option>
                        <option  value="500">500</option>
                        <option  value="1000">1000</option>
                        <option  value="10000">10000</option>
                        <option selected="selected"  value="50000">50000</option>
                        <option  value="100000">100000</option>
                      </select>                       
                    </span>
                    <a class="helplink hiding" title="help" id="hitsizeHelp" href="#"><i class="fas fa-question-circle"></i> <span class="usa-sr-only">Help</span></a>
                    <p toggle="hitsizeHelp" class="helpbox hidden">
                      Maximum number of database sequences (with unique sequence identifier) Blast finds for primer-blast to screen for primer pair specificities. Note that the actual number of similarity regions (or the number of hits) may be much larger than this (for example, there may be a large number of hits on a single target sequence such as a chromosome).   Choose a higher value if you need to perform more stringent search.
                    </p>      

```
You can find the description and suggested values for this search option. HITSIZE ='1000000' is added to the search below along with several options that increase the number of entries returned from primer_search.


### Example
```
forward_primer_seq = "TAGAACAGGCTCCTCTAG"
reverse_primer_seq =  "TTAGATACCCCACTATGC"
output_directory_path <- "/my/directory/12S_V5F1_remote_111122_modified_params"
metabarcode_name <- "12S_V5F1"
accession_taxa_sql_path <- "/my/directory/accessionTaxa.sql"
blast_db_path <- "/my/directory/ncbi_nt/nt"

get_seeds_remote(forward_primer_seq,
                reverse_primer_seq,
                output_directory_path,
                metabarcode_name,
                accession_taxa_sql_path,
                HITSIZE ='1000000',
                evalue='100000',
                word_size='6',
                MAX_TARGET_PER_TEMPLATE = '5',
                NUM_TARGETS_WITH_PRIMERS ='500000', minimum_length = 50,
                MAX_TARGET_SIZE = 200,
                organism = c("1476529", "7776"), return_table = FALSE)

# This results in approximately 111500 blast seed returns (there is some variation due to database updates, etc.), note the default generated approximately 1047.
# This assumes the user is not throttled by memory limitations.              
```


## [blast_seeds](https://limey-bean.github.io/blast_seeds)

<img src="/flowcharts/blast_seeds-flowchart.png" width = 10000 />

### Overview
blast_seeds uses [blast_datatable](https://limey-bean.github.io/blast_datable) to search against a blast
formatted database. It creates a permanent directory `blast_seeds_output` and
a temporary directory `blast_seeds_save` in the `output_directory_path`. It saves from
and passes files to [blast_datatable](https://limey-bean.github.io/blast_datable) while the run is in progress.
During the final steps of the function the final data is saved in
`rblast_seeds_output` recording the results of the blast.

The [blast_datatable](https://limey-bean.github.io/blast_datable) call saves intermediate results and
metadata about the search as local files in the save directory generated by
blast_seeds. This allows the function to resume a partially
completed blast, mitigating the consequences of encountering an
error or experiencing other interruptions. To resume a partially completed
blast, supply the same seeds and working directory. See the documentation
of [blast_datatable](https://limey-bean.github.io/blast_datable) for more information.

### Expected Output
During the blast_seeds the following data are cashed as files and passed into
[blast_datatable](https://limey-bean.github.io/blast_datable): output_table.txt (most recent updates from the
blast run), blast_seeds_passed_filter.txt (seed table that tracks the blast
status of seeds), unsampled_indices.txt (list of seed indices that need to
be blasted), too_many_ns.txt (tracks seeds that have been removed due to more
consecutive Ns in a sequence than are acceptable (see parameter wildcards),
blastdbcmd_failed.txt (tracks reads that are present in the seeds database,
but not the local blast database. This is relevant for the results of
[get_seeds_remote](https://limey-bean.github.io/get_seeds_remote), and lastly num_rounds.txt (tracks the number of
completed blast round for a given seed file).

The final output of blast_seeds, stored in `blast_seeds_output`, are the
following: summary.csv (blast output with appended taxonomy),
{metabarcode_name}_.fasta, {metabarcode_name}.taxonomy,
{metabarcode_name}_blast_seeds_summary_unique_taxonomic_rank_counts.txt,
too_many_ns.txt, blastdbcmd_failed.txt.

### Detailed Steps

[blast_seeds](https://limey-bean.github.io/blast_seeds) passes a datatable returned by [get_seeds_remote](https://limey-bean.github.io/get_seeds_remote), or [get_seeds_local](https://limey-bean.github.io/get_seeds_local), [blast_datatable](https://limey-bean.github.io/blast_datable) which uses a random stratified sample based on taxonomic rank to iteratively blast and process the seeds in the datatable. The user can specify how many sequences can be blasted simultaneously using max_to_blast. The randomly sampled seeds (or subsets of seeds) are sent to [run_blastdbcmd_blastn_and_aggregate_resuts](https://limey-bean.github.io/run_blastdbcmd_blastn_and_aggregate_resuts), which uses [run_blastdbcmd](https://limey-bean.github.io/run_blastdbcmd)to find a seed sequence that corresponds to the accession number and forward and reverse stops recorded in the seeds table. [run_blastdbcmd](https://limey-bean.github.io/run_blastdbcmd) outputs sequences as .fasta-formatted strings, which
[run_blastdbcmd_blastn_and_aggregate_resuts](https://limey-bean.github.io/run_blastdbcmd_blastn_and_aggregate_resuts) concatenates into a multi-line fasta, then passes to [run_blastn](https://limey-bean.github.io/run_blastn) as an argument. The output of
[run_blastn](https://limey-bean.github.io/run_blastn) is de-replicated by accession, and only the longest read per replicates is retained in the output table. The run state is saved and passed back to [blast_datatable](https://limey-bean.github.io/blast_datable).

Once all of the seeds of the random sample are processed, they are removed from the dataframe as are the seeds found as blast hits. blast-datatable repeats this process or stratified random sampling until there are fewer reads remaining than max_to_blast, at which point it blasts all remaining seeds. The final aggregated results are cleaned for multiple blast taxids, hyphens, and wildcards and returned with taxonomy added using [get_taxonomizr_from_accession](https://limey-bean.github.io/get_taxonomizr_from_accession). Additionally, it saves its state as text files in a specified directory with each iteration, allowing the user to restart an interrupted run of [blast_seeds](https://limey-bean.github.io/blast_seeds).

**Note:**
The blast db downloaded from NCBIs FTP site has representative
accessions meaning identical sequences have been collapsed across multiple
accessions even if they have different taxid. Here we identify representative
accessions with multiple taxids, and unpack all of the accessions that were
collapsed into that representative accessions. [blast_seeds](https://limey-bean.github.io/blast_seeds) does not identify or unpack representative accessions that report a single taxid.

**Saving data:**
blast_datatable uses files generated in [run_blastdbcmd_blastn_and_aggregate_resuts](https://limey-bean.github.io/run_blastdbcmd_blastn_and_aggregate_resuts) that store intermediate
results and metadata about the search to local files as it goes. This allows
the function to resume a partially completed blast, partially mitigating
the consequences of encountering an error or experiencing other interruptions.
Interruptions while blasting a subset of a random stratified sample will
result in a loss of the remaining reads of the subsample, and may decrease
overall blast returns. The local files are written to `blast_seeds_save` by
[rsave_state](https://limey-bean.github.io/save_state). Manually changing these files is not suggested as
it can change the behavior of blast_datatable.

**Restarting an interrupted [blast_seeds](https://limey-bean.github.io/blast_seeds) run:**
To restart from an incomplete blast_datatable, submit the previous command
again. Do not modify the paths specified in the previous command, however
parameter arguments (e.g. rank, max_to_blast) can be modified. blast_datable
will automatically detect save files and resume from where it left off.

**Warning:**
If you are resuming from an interrupted blast, make sure you supply
the same data.frame for [blast_seeds](https://limey-bean.github.io/blast_seeds). If you intend to start a new blast,
make sure that there is not existing blast save data in the directory you
supply for `blast_seeds_save`.

**Note:**
[blast_datatable](https://limey-bean.github.io/blast_datable) does not save intermediate data
from [run_blastdbcmd](https://limey-bean.github.io/run_blastdbcmd), so if it is interrupted while getting building the fasta to
submit to [run_blastn](https://limey-bean.github.io/run_blastn) it will need to repeat some work when resumed. The argument
`max_to_blast` controls the frequency with which it calls blastn, so it can
be used to make [blast_datatable](https://limey-bean.github.io/blast_datable) save more frequently.

### Parameters
**seeds_output_path**
+ a path to a csv from get_seeds_local or get_seeds_remote
+         e.g. seeds_output_path <- '/my/rCRUX_output_directory/12S_V5F1_filtered_get_seeds_remote_output_with_taxonomy.csv'
**blast_db_path**
+ a directory containing a blast-formatted database
+         e.g blast_db_path <- "/my/ncbi_nt/nt"
**accession_taxa_sql_path**
+ a path to the accessionTaxa sql created by taxonomizr.
+         e.g. accession_taxa_sql_path <- "/my/accessionTaxa.sql"
**output_directory_path**
+ a directory in which to save partial and complete output.
+         e.g. output_directory_path = "/path/to/output/12S_V5F1_local_111122_e300_111122"
**metabarcode_name**
+ a prefix for the output fasta, taxonomy, and count of unique ranks.
+         e.g. metabarcode_name <- "12S_V5F1"
**expand_vectors**
+ logical, determines whether to expand too_many_Ns
        and not_in db into real tables and write them in the output directory.
+        The default is expand_vectors = TRUE
**warnings**
+ value to set the "warn" option to during the function call.
        On exit it returns to the previous value. Setting this argument to
        NULL will not change the option.
**...**
+ additional arguments passed to [blast_datatable](https://limey-bean.github.io/blast_datatable)
**sample_size**
+ passed to [blast_datatable](https://limey-bean.github.io/blast_datatable) is the the number of
        entries to sample per rank.
+       The default sample_size = 1 - is recommended
+       Note: unless the user is sampling higher order taxonomy.  If there are not enough seeds to sample per rank the run will end in an error.
**max_to_blast**
+ passed to [blast_datatable](https://limey-bean.github.io/blast_datatable) and is the maximum
        number of entries to accumulate into a fasta before calling blastn.
+        The default is max_to_blast = 1000
+        Note: the optimal number of reads to blast will depend on the user's environment (available RAM) and the number of possible hits (determined by marker and parameters)
**wildcards**
+ passed to [blast_datatable](https://limey-bean.github.io/blast_datatable) us a character vector
        that represents the minimum number of consecutive Ns the user will
        tolerate in a given seed or hit sequence.
+       The default iswildcards = "NNNN"
**rank**
+ passed to [blast_datatable](https://limey-bean.github.io/blast_datatable) is the data column
        representing the taxonomic rank to randomly sample.
+       The default is rank = genus
+       Note: sampling a lower rank  (e.g. species) will generate more total hits and take more time, conversely sampling a higher rank (e.g. family) will generate fewer total hits and take less time.
**ncbi_bin**
+ passed to [run_blastdbcmd](https://limey-bean.github.io/run_blastdbcmd)
          and [run_blastn](https://limey-bean.github.io/run_blastn)is
        the path to blast+ tools if not in the user's path.  Specify only if
        blastn and blastdbcmd  are not in your path.
+         The default is ncbi_bin = NULL
+         Note: if not specified in path do the following: ncbi_bin = "/my/local/blast+_folder/bin".
**evalue**
+ passed to [run_blastn](https://limey-bean.github.io/run_blastn) is the number of expected hits
       with a similar quality score found by chance.
+      The default is evalue = 1e-6
**coverage****
+ passed to [run_blastn](https://limey-bean.github.io/run_blastn) is the minimum percent of the
        query length recovered in the subject hits.
+       The default is coverage = 50
**perID**
+ passed to [run_blastn](https://limey-bean.github.io/run_blastn) is the minimum percent identity
        of the query relative to the subject hits.
+       The default is perID = 70
**align**
+ passed to [run_blastn](https://limey-bean.github.io/run_blastn) is the maximum number of subject
        hits to return per query blasted.
+        The default is align = 50000


### Example
```
seeds_output_path <- "/my/directory/12S_V5F1_remote_111122_modified_params/blast_seeds_output/summary.csv""
output_directory_path <- "/my/directory/12S_V5F1_remote_111122_modified_params"
metabarcode_name <- "12S_V5F1"
accession_taxa_sql_path <- "/my/directory/accessionTaxa.sql"
blast_db_path <- "/my/directory/ncbi_nt/nt"


blast_seeds(seeds_output_path,
            blast_db_path,
            accession_taxa_sql_path,
            output_directory_path,
            metabarcode_name,
            rank = species,
            max_to_blast = 750)

# using the rank of species will increase the number of total unique blast hits
# modifying the max_to_blast submits fewer reads simultaneously and reduces overall RAM while extending the run

```


## [derep_and_clean_db](https://limey-bean.github.io/derep_and_clean_db)

<img src="/flowcharts/derep_and_clean_db-flowchart.png" width = 10000 />

### Overview
derep_and_clean_db takes the output from [blast_seeds](https://limey-bean.github.io/blast_seeds) and de-replicates the dataset to identify representative sequences.

### Expected Outputs
It generates an output directory called `derep_and_clean_db` at `output_directory_path` to store the output .csv files and the fasta and taxonomy file generated by the function.


### Detailed Steps
Before de-replicating the data set, all sequences with NA taxonomy for phylum,
class, order, family, and genus are removed from the dataset because they
typically represent environmental samples with low value for taxonomic classification. These sequences are stored in a
`Sequences_with_mostly_NA_taxonomic_paths.csv`

All sequences with the same length and composition are collapsed to a single
database entry, where the accessions and taxids (if there are more than one)
are concatenated. The sequences with a clean taxonomic path (e.g. no ranks
with multiple entries) are written to `Sequences_with_single_taxonomic_path.csv`.

Sequences with multiple entries for a given taxonomic rank are written to
`Sequences_with_multiple_taxonomic_paths.cvs`. These sequences are processed
further by removing NAs from rank instances with more than one entry
(e.g. "Chordata, NA" will mutate to "Chordata"). Any remaining instances
of taxonomic ranks more than one taxid will be reduced to NA
(e.g. "Badis assamensis, Badis badis" will mutate to "NA"). These sequences,
with taxonomic paths shortened to the lowest taxonomic agreement, are written
to `Sequences_with_lowest_common_taxonomic_path_agreement.csv`.

Lastly, the sequences from `Sequences_with_single_taxonomic_path.csv`
and `Sequences_with_lowest_common_taxonomic_path_agreement.csv` are used to
generate a fasta file and taxonomy file of representative NCBI accessions for
each sequence.  The number of accessions identical to the representative
accession is given.

### Parameters
**output_directory_path**
+ the path to the output directory
+        e.g. "/path/to/output/12S_V5F1_remote_111122"

**summary_path**
+ the path to the input file
+        e.g. "/path/to/output/12S_V5F1_remote_111122/blast_seeds_output/summary.csv"
**metabarcode_name**
+ used to name the subdirectory and the files.
+        e.g. metabarcode_name <- "12S_V5F1"

### Example
```
output_directory_path <- "/my/directory/12S_V5F1_remote_111122_modified_params"
summary_path <- "/my/directory/12S_V5F1_remote_111122_modified_params/blast_seeds_output/summary.csv"
metabarcode_name <- "12S_V5F1"


derep_and_clean_db(output_directory_path, summary_path, metabarcode_name)

```

## References

<div id="refs" class="references csl-bib-body hanging-indent">

</div>

Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.

</div>

Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) "BLAST+: architecture and applications." BMC Bioinformatics 10:421.

</div>


<div id="ref-Curdetal.2019" class="csl-entry">

Curd, E.E., Gold, Z., Kandlikar, G.S., Gomer, J., Ogden, M., O'Connell, T.,
Pipes, L., Schweizer, T.M., Rabichow, L., Lin, M. and Shi, B., 2019.
Anacapa Toolkit: An environmental DNA toolkit for processing multilocus
metabarcode datasets. Methods in Ecology and Evolution, 10(9), pp.1469-1475.
<https://doi.org/10.1111/2041-210X.13214>.

</div>

Hester, J., 2020. primerTree: Visually Assessing the Specificity and Informativeness of Primer Pairs.

</div>
<div id="ref-R Core Team 2018" class="csl-entry">
R Core Team, R., 2021. R: A language and environment for statistical computing.

</div>
</div>

Sherrill-Mix, S., 2019. taxonomizr: Functions to Work with NCBI Accessions and Taxonomy. R package version 0.5. 3.

</div>
<div id="ref-Ye et al. 2012" class="csl-entry">

Ye J, Coulouris G, Zaretskaya I, Cutcutache I, Rozen S, & Madden TL. (2012) "Primer-BLAST: a tool to design target-specific primers for polymerase chain reaction." BMC Bioinformatics 13:134.

</div>
