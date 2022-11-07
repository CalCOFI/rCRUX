# rCRUX: Generate CRUX metabarcoding reference libraries in R



**Authors:** [Luna Gal](<https://github.com/LunaGal>), [Zachary Gold](<https://github.com/zjgold>), [Ramon Gallego](https://github.com/ramongallego), [Emily Curd](<https://github.com/limey-bean>)<br/>
**Inspiration:**
The late, great [Jesse Gomer](<https://github.com/jessegomer?tab=repositories>). Coding extraordinaire and dear friend.<br/>
**License:**
[GPL-3](https://opensource.org/licenses/GPL-3.0) <br/>
**Support:**
Support for the development of this tool was provided by [CalCOFI](<https://calcofi.org/>) and NOAA. <br/>

eDNA metabarcoding is increasingly used to survey biological communities using common universal and novel genetic loci. There is a need for an easy to implement computational tool that can generate metabarcoding reference libraries for any locus, and are specific, and comprehensive. We have reimagined CRUX [Curd et al. 2019](https://doi.org/10.1111/2041-210X.13214) and developed the rCRUX package to fit this need by generating taxonomy and fasta files for any user defined locus.  The typical workflow involves using get_seeds_local() or get_seeds_remote() to download and wrangle results from NCBI’s primer BLAST tool, then using blast_seeds() to search a local NCBI database for matches.


## Typical Workflow
<img src="RCRUX_Flowchart.pdf" width = 10000 />


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

#### NOTE: These only need to be downloaded once or as NCBI updates databases. </br>

**BLAST+**

NCBI's [BLAST+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/) suite must be locally installed and accessible in the user's path. NCBI provides installation instructions for [Windows](https://www.ncbi.nlm.nih.gov/books/NBK52637/), [Linux](https://www.ncbi.nlm.nih.gov/books/NBK52640/), and [Mac OS](https://www.ncbi.nlm.nih.gov/books/NBK569861/). Version 2.10.1+ is verified compatible with rCRUX.

**Blast-formatted database**

rCRUX requires a local blast-formatted nucleotide database. These can be user generated or download a pre-formatted database from [NCBI](https://ftp.ncbi.nlm.nih.gov/blast/db/).  NCBI provides a tool (perl script) for downloading databases as part of the blast+ package. A brief help page can be found [here](https://www.ncbi.nlm.nih.gov/books/NBK569850/).

The following shell script can be used to download the blast-formatted nucleotide database.

```
mkdir NCBI_blast_nt
cd NCBI_blast_nt
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt*
time for file in *.tar.gz; do tar -zxvf $file; done
cd ..

```

You can test your nt blast database using the following command:

```
blastdbcmd -db '/Users/limeybean/Dropbox/CRUX_2.0/ncbi_nt/nt' -dbtype nucl -entry MN937193.1 -range 499-633
```
If you do not get the following, something went wrong in the build.

```
>MN937193.1:499-633 Jaydia carinatus mitochondrion, complete genome
TTAGATACCCCACTATGCCTAGTCTTAAACCTAGATAGAACCCTACCTATTCTATCCGCCCGGGTACTACGAGCACCAGC
TTAAAACCCAAAGGACTTGGCGGCGCTTCACACCCACCTAGAGGAGCCTGTTCTA
```

Possible error include but are not limited to:
1. Partial downloads of database files. Extracting each TAR archive (e.g. nt.00.tar.gz.md5) should cresult in 8 files with the following extensions(.nhd, .nhi, .nhr, .nin, .nnd, .nni, .nog, and .nsq).  If a few archives fail during download, you can redownload and unpack those that failed. You do not have to redownload all archives.  

2. You downloaded and built a blast database from ncbi fasta files but did not specify -parse_seqids


The nt database is **~242 GB** (as of 8/31/22) and can take several hours (overnight) to build. Loss of internet connection can lead to partially downloaded files and blastn errors (see above). rCRUX can access and successfully build metabarcode references using databases stored on external drives.

**Taxonomizr**

rCRUX uses the [taxonomizr](https://cran.r-project.org/web/packages/taxonomizr/vignettes/usage.html) package for taxonomic assignment based on NCBI [Taxonomy id's \(taxids\)](https://www.ncbi.nlm.nih.gov/). Many rCRUX functions require a path to a local taxonomizr readable sqlite database. This database can be built using taxonomizr's [prepareDatabase](https://www.rdocumentation.org/packages/taxonomizr/versions/0.8.0/topics/prepareDatabase) function.

This database is **~72 GB** (as of 8/31/22) and can take several hours (overnight) to build. Loss of internet connection can lead to partially downloaded files and taxonomizr run errors. rCRUX can access and successfully build metabarcode references using databases stored on external drives.

The following code can be used to build this database:

```
library(taxonomizr)

accession_taxa_sql_path <- "/my/accessionTaxa.sql"
prepareDatabase(accession_taxa_sql_path)

```

**Note** Please see the taxononmizr readme for manual installation of the accessionTaxa.sql database. This can be required with poor bandwidth connections. If built manually, make sure to delete any files other than the accessionTaxa.sql database (e.g. keeping nucl_gb.accession2taxid.gz leads to a warning message).

# Example pipeline

The following example shows a simple rCRUX pipeline from start to finish. Note that this example will require internet access and considerable database storage (~314 GB, see section above), run time (mainly for blastn), and system resources to execute.

* Note that local blast and taxonomic assignment databases can be stored on an external hard drive. It increases run time, but is a good option if computer storage capacity is limited.

There are two options to generate seeds for the database generating blast step [rCRUX::blast_seeds()]: get_seeds_local and get_seeds_remote. The local option is slower, however it is not subject to the memory limitations of using the NCBI primer_blast API. The local option is recommended if the user is building large database, wants to include any taxid in the search, and has many degenerate sites in their primer set.

**get_seeds_local**


This example uses default parameters, with the exception of evalue, to minimize run time.

```
output_directory_path <- "/my/rCRUX_output_directory/12S_V5F1_default_params" #path to desired output directory

metabarcode_name <- "12S_V5F1" # desired name of metabarcode locus

accession_taxa_sql_path <- "/my/accessionTaxa.sql"

forward_primer_seq = "TAGAACAGGCTCCTCTAG"

reverse_primer_seq =  "TTAGATACCCCACTATGC"

blast_db_path <- "/my/ncbi_nt/nt"


get_seeds_local(forward_primer_seq,
                 reverse_primer_seq,
                 output_directory_path,
                 metabarcode_name,
                 accession_taxa_sql_path,
                 blast_db_path, evalue = 300)

```

Two output .csv files are automatically created at this path based on the arguments passed to get_seeds_local.  One includes taxonomy the other does not.

A unique taxonomic rank summary file is also generated (e.g. the number of unique phyla, class, etc in the blast hits). If a taxonomic rank category contains NA's, they will be counted as a single unique rank. Sequence availability in NCBI for a given taxid is a limiting factor.


**get_seeds_remote**

This example uses default parameters to minimize run time.

Searching jawless vertebrates (taxid: "1476529") and jawed vertebrates (taxid: "7776").

```
output_directory_path <- "/my/rCRUX_output_directory/12S_V5F1_default_params" #path to desired output directory

metabarcode_name <- "12S_V5F1" # desired name of metabarcode locus

accession_taxa_sql_path <- "/my/accessionTaxa.sql"

forward_primer_seq = "TAGAACAGGCTCCTCTAG"

reverse_primer_seq =  "TTAGATACCCCACTATGC"


get_seeds_remote(forward_primer_seq,
          reverse_primer_seq,
          output_directory_path,
          metabarcode_name,
          accession_taxa_sql_path,
          organism = c("1476529", "7776"),
          return_table = FALSE)

```

*Note that using default parameters only 1047 hits are returned from NCBI's primer blast (as of 9-25-22).*

Two output .csv files are automatically created at this path based on the arguments passed to get_seeds_remote.  One includes taxonomy the other does not.

A unique taxonomic rank summary file is also generated (e.g. the number of unique phyla, class, etc in the blast hits). If a taxonomic rank category contains NA's, they will be counted as a single unique rank.

Sequence availability in NCBI for a given taxid is a limiting factor.

[Modifying defaults can increase the number of returns by orders of magnitude.](#Search-options)




**blast_seeds**


Iterative searches are based on randomly sampling unique taxonomic groups for a given rank from the get_seeds_local output table to create a set of blast seeds. For example, the default is to randomly sample one read from each genus.  The user can select any taxonomic rank present in the get_seeds_local output table. The number of seeds selected may exceed the users available RAM, and for that reason the user can choose the maximum number of reads to blast at one time (max_to_blast, default = 1000). blast_seeds will subsample each set of seeds based on max_to_blast and process all seeds before starting a new search for seeds to blast.


```
seeds_output_path <- '/my/rCRUX_output_directory/12S_V5F1_primerTree_output_with_taxonomy.csv'
# this is output from get_seeds_local or get_seeds_remote

blast_db_path <- "/my/local/blast_database/nt"

accession_taxa_sql_path <- "/my/accessionTaxa.sql"

output_directory_path <- '/my/rCRUX_output_directory/12S_V5F1_default_params/'

metabarcode_name <- "12S_V5F1"


blast_seeds(seeds_output_path,
            blast_db_path,
            accession_taxa_sql_path,
            output_directory_path,
            metabarcode_name)    
```


After each round of blast, the system state is saved. If the script is terminated after a round of blast, the user can pick up where they left off. The user can also change parameters at this point (e.g. change the max_to_blast or rank)

The output includes a summary table of unique blast hits, a multi fasta file, a taxonomy file, a unique taxonomic rank summary file, a list of all of the accessions not present in your blast database, and a list of accessions with 4 or more Ns in a row (default for that parameter is wildcards = "NNNN") the default number of reads to blast per rank is 1 (default for that parameter is sample_size = 1). The script will error out if the user asks for more reads per rank than exist in the blast seeds table.     



**If BLAST+ is not in your path do the following**


```
blast_seeds(seeds_output_path,
            blast_db_path,
            accession_taxa_sql_path,
            output_directory_path,
            metabarcode_name,
            ncbi_bin = "/my/local/blast+_folder")


```

Example output can be found [here](/examples/12S_V5F1_generated_9-21-22).

**Note**, there will be variability between runs due to primer blast return parameters and random sampling of the blast seeds table that occurs during blast_seeds.

# Detailed Explanation of The Major Functions

## Internal data pipeline

# [get_seeds_local](https://lunagal.github.io/get_seeds_local)

This script takes a set of forward and reverse primer sequences and generates csv summaries of data returned from a locally run adaptation of [NCBI's primer blast](https://www.ncbi.nlm.nih.gov/tools/primer-blast/) to find possible full length barcode sequence containing forward and reverse primer matches. It also generates a count of unique instances of taxonomic ranks (Phylum, Class, Order, Family, Genus, and Species)

add info....

## Search options
add info...


# [get_seeds_remote](https://lunagal.github.io/get_blast_seeds)

This script takes a set of forward and reverse primer sequences and generates csv summaries of data returned from [NCBI's primer blast](https://www.ncbi.nlm.nih.gov/tools/primer-blast/) about full length barcode sequence containing primer matches. It also generates a count of unique instances of taxonomic ranks (Phylum, Class, Order, Family, Genus, and Species)

get_seeds_remote uses modified versions of functions from the [primerTree](https://CRAN.R-project.org/package=primerTree) package to submit queries to NCBI's primer BLAST tool, then aggregates results into a single data.frame. primer_search (Modified from [primerTree](https://CRAN.R-project.org/package=primerTree)) expands degenerate primers into each possible non-degenerate primer and submits a query for each. get_seeds_remote further multiplies the number of queries by allowing the user to query the primers for each organism in a vector. get_seeds_remote collects all these results from primer_search, filters them based on product length, and adds taxonomic data using the taxonomizr package.

### Organism(s)
primer BLAST defaults to homo sapiens, so it is important that you supply a specific organism or organisms. NCBI's taxids can be found [here](https://www.ncbi.nlm.nih.gov/taxonomy). You can specify multiple organism by passing a character vector containing each of the options, like in the example below.

## Search options

get_seeds_remote passes many parameters to NCBI's primer blast tool. You can match the parameters to the fields available in the GUI here. First, use your browser to view the page source. Search for the field you are interested in by searching for the title of the field. It should be enclosed in a tag. Inside the label tag, it says for = "<name_of_parameter>". Copy the string after for = and add it to get_seeds_remote as the name of a parameter, setting it equal to whatever you like.

As of 2022-08-16, the primer blast GUI contains some options that are not implemented by primer_search. The [table below](#Table-of-available-options) documents available options.

**You can checking [primerblast](https://www.ncbi.nlm.nih.gov/tools/primer-blast/) for more information on how to modify search options.  For example, if you to generate a larger hitsize, open the source of the primer designing tool and look for that string. You find the following:

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

```
output_directory_path <- "/Users/limeybean/Dropbox/CRUX_2.0/12S_V5F1_modified_params"

metabarcode_name <- "12S_V5F1"

accession_taxa_sql_path <- "/Users/limeybean/Dropbox/CRUX_2.0/accession2taxid/accessionTaxa.sql"


forward_primer_seq = "TAGAACAGGCTCCTCTAG"

reverse_primer_seq =  "TTAGATACCCCACTATGC"

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
```

This results in 111500 blast seed returns, note the default generated 1047. This assumes the user is not throttled by memory limitations.


Example output can be found [here](/examples/12S_V5F1_generated_9-21-22).

### Table of available options

**Need to check for accuracy and completeness**

| Name                                   | rCrux Default  |
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




# [blast_seeds](https://lunagal.github.io/blast_seeds)

blast_seeds uses the entries generated by get_seeds_remote or get_seeds_local and the nucleotide-nucleotide matching of blastn to generate a .csv of NCBI database entries that match a sequence found in the get_seeds_local or get_seeds_remote step.

blastseeds is a wrapper function that passes data to blast_datatable. blast_seeds handles the creation of a hidden save directory and an output directory and writes a .csv summarizing the results of blast_datatable.  It also generates a multi fasta file and corresponding a taxonomy file for all recovered hits. A unique taxonomic rank summary file (same format as the get blast seeds output file), a list of all of the accessions not present in your blast database, and a list of accessions with 4 or more Ns in a row.

Internally, blast_datatable iteratively samples rows from the table of blast seeds, calls blastdbcmd on each accession number, and uses blastn to build a table of nucleotide matches. It samples by drawing random indices from a list of unsampled indices and examining the rows at those indices. Specifically, it randomly samples single indices from unique instances of a taxonomic rank (e.g. one accession number per genus). Depending on the number of seeds to be processed in each iteration, the next step may be run on subsamples of the seed set.

blast_datatable passes those seeds to run_blastdbcmd, which extracts the accession number, forward stop, and reverse stop, then uses them as arguments for blastdbcmd. blastdbcmd outputs a fasta, which blast_datatable aggregates into a multi-fasta character vector. blast_datatable purges any entry that has more than a specified number of Ns or did not return a result, recording those indices. When it has finished building the mutli-fasta vector, it passes it to blastn, which returns every nucleotide sequence that matches a sequence in the file. run_blastn parses the blastn output into a data.frame, and blast_datatable adds that data.frame to its output. All output is deduplicated by accession retaining the duplicate with the longest sequence.

After each iteration, the accessions recovered through blastn are removed from the table of blast seeds. Then the process is repeated until every row is sampled. However, once the number of seeds in the table reaches the max_to_blast parameter, all remaining seeds will be treated as a set of blast seeds.  

After all blast seeds are processed, taxonomizr is used to add taxonomic data to the data.frame based on the accession numbers. The final output is the aggregate of all blastn calls with the taxonomic data added.


## Funding

The [CalCOFI](https://calcofi.org/) program provided funding support for this project.

## Acknowledgments

This work benefited from the amazing input of many including Lenore Pipes, Sarah Stinson, Gaurav Kandlikar, and Maura Palacios Mejia.

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Curdetal.2019" class="csl-entry">

Curd, E.E., Gold, Z., Kandlikar, G.S., Gomer, J., Ogden, M., O'Connell, T.,
Pipes, L., Schweizer, T.M., Rabichow, L., Lin, M. and Shi, B., 2019.
Anacapa Toolkit: An environmental DNA toolkit for processing multilocus
metabarcode datasets. Methods in Ecology and Evolution, 10(9), pp.1469-1475.
<https://doi.org/10.1111/2041-210X.13214>.

</div>

<div id="ref-next" class="csl-entry">

next

</div>

</div>
