# rCRUX: Generate CRUX metabarcoding reference libraries in R



**Authors:** [Luna Gal](<https://github.com/LunaGal>), [Zachary Gold](<https://github.com/zjgold>), [Emily Curd](<https://github.com/limey-bean>)<br/>
**License:**
[GPL-3](https://opensource.org/licenses/GPL-3.0)


eDNA metabarcoding is increasingly used to survey biological communities using common universal and novel genetic loci. There is a need for an easy to implement computational tool that can generate metabarcoding reference libraries for any locus, and are specific, and comprehensive. We have reimagined CRUX [Curd et al. 2019](https://doi.org/10.1111/2041-210X.13214) and developed the rCRUX package to fit this need by generating taxonomy and miltifasta files for any user defined locus.  The typical workflow involves using get_blast_seeds() to download and wrangle results from NCBI’s primer BLAST tool, then using rcrux_blast() to search a local NCBI database for matches.


## Typical Workflow
Needs to be updated, but general flow is still accurate.
<img src="RCRUX_Flowchart.png" width = 10000 />


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

**BLAST+**

NCBI's [BLAST+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/) suite must be locally installed and accessible in the user's path. NCBI provides installation instructions for [Windows](https://www.ncbi.nlm.nih.gov/books/NBK52637/), [Linux](https://www.ncbi.nlm.nih.gov/books/NBK52640/), and [Mac OS](https://www.ncbi.nlm.nih.gov/books/NBK569861/). Version 2.10.1+ is verified compatible with rCRUX.

**Blast-formatted database**

rCRUX requires a local blast-formatted nucleotide database. These can be user generated or download a pre-formatted database from [NCBI](https://ftp.ncbi.nlm.nih.gov/blast/db/).  NCBI provides a tool (perl script) for downloading databases as part of the blast+ package. A brief help page can be found [here](https://www.ncbi.nlm.nih.gov/books/NBK569850/).

The nt database is **~242 GB** (as of 8/31/22) and can take several hours (overnight) to build. Loss of internet connection can lead to partially downloaded files and blastn errors. rCRUX can access and successfully build metabarcode references using databases stored on external drives.

**Taxonomizr**

rCRUX uses the [taxonomizr](https://cran.r-project.org/web/packages/taxonomizr/vignettes/usage.html) package for taxonomic assignment based on NCBI [Taxonomy id's \(taxids\)](https://www.ncbi.nlm.nih.gov/). Many rCRUX functions require a path to a local taxonomizr readable sqlite database. This database can be built using taxonomizr's [prepareDatabase](https://www.rdocumentation.org/packages/taxonomizr/versions/0.8.0/topics/prepareDatabase) function.

This database is **~72 GB** (as of 8/31/22) aand can take several hours (overnight) to build. Loss of internet connection can lead to partially downloaded files and taxonomizr run errors. rCRUX can access and successfully build metabarcode references using databases stored on external drives.

The following code can be used to build this database:

```
library(taxonomizr)

accession_taxa_path <- "/my/accessionTaxa.sql"
prepareDatabase(accession_taxa_path)

```


# Example pipeline

The following example shows a simple rCRUX pipeline from start to finish. Note that this example will require internet access and considerable database storage (see section above), run time (mainly for blastn), and system resources to execute.

* Note that local blast and taxonomic assignment can be run on databases stored on an external hard drive. It increases run time, but is a good option if computer storage capacity is limited.

**get_blast_seeds**
Searching jawless vertebrates (taxid: "1476529") and jawed vertebrates (taxid: "7776").

This example uses default parameters to minimize run time.
```
blast_seeds_parent <- "/my/rCRUX_output_directory/12S_V5F1_default_params"
metabarcode <- "12S_V5F1"
accession_taxa_path <- "/my/accessionTaxa.sql"


get_blast_seeds("TAGAACAGGCTCCTCTAG", "TTAGATACCCCACTATGC",
                 blast_seeds_parent, metabarcode, accession_taxa_path,
                 organism = c("1476529", "7776"), return_table = FALSE)


# Output .csv files are automatically created at this path based on the arguments passed to get_blast_seeds
# A unique taxonomic rank summary file is also generated. If a taxonomic rank catergory contains NA's, they will be counted as a single unique rank.
# Note that using default parameters only 1047 hits are returned from NCBI's primer blast.  
# Sequence availability in NCBI for a given taxid is a limiting factor.
```
[Modifying defaults can increase the number of returns by orders of magnitude.](#Search-options)

**rcrux_blast**
Searches are based on randomly sampling unique taxonomic groups for a given rank from the get_blast_seeds output table. For example, the default is to randomly sample one read from each genus.  The user can select any taxonomic rank present in the get_blast_seeds output table, and should choose a rank that can be feasibly run in their environment (e.g. 20K + reads are too memory intensive for many laptops, and some targets return large numbers of hits that can max out RAM), and provide the user sufficient resolution. If there are 1,000 (or user defined max_to_blast) or fewer blast seeds to process, for a given round of blasting, the entire blast seeds table will be blasted.
```
seeds_path <- '/my/rCRUX_output_directory/12S_V5F1_primerTree_output_with_taxonomy.csv'
# this is output from get_blast_seeds
db_dir <- "/my/local/blast_database/nt"
accession_taxa_path <- "/my/accessionTaxa.sql"
working_dir <- '/my/rCRUX_output_directory/12S_V5F1_default_params/'
metabarcode <- "12S_V5F1"


rcrux_blast(seeds_path, db_dir, accession_taxa_path, working_dir,
            metabarcode)

# the default number of reads to blast per rank is 1. The script will error out if the user asks for more reads per rank than exist in the blast seeds table.         
```
**If BLAST+ is not in your path do the following**
```
rcrux_blast(seeds_path, db_dir, accession_taxa_path, working_dir,
            metabarcode, ncbi_bin = "/my/local/blast+_folder")


```

Example output can be found [here](/examples/12S_V5F1_generated_9-21-22).

**Note**, there will be variability between runs due to primer blast return parameters and random sampling of the blast seeds table that occurs during rcrux_blast.

# Detailed Explaination of The Major Functions

# [get_blast_seeds](https://lunagal.github.io/get_blast_seeds)

This script takes a set of forward and reverse primer sequences and generates csv summaries of data returned from [NCBI's primer blast](https://www.ncbi.nlm.nih.gov/tools/primer-blast/) about full length barcode sequence containing primer matches. It also generates a count of unique instances of taxonomic ranks (Phylum, Class, Order, Family, Genus, and Species)

get_blast_seeds uses modified versions of functions from the [primerTree](https://CRAN.R-project.org/package=primerTree) package to submit queries to NCBI's primer BLAST tool, then aggregates results into a single data.frame. primer_search expands degenerate primers into each possible non-degenerate primer and submits a query for each. get_blast_seeds further multiplies the number of queries by allowing the user to query the primers for each organism in a vector. get_blast_seeds collects all these results from primer_search, filters them based on product length, and adds taxonomic data using the taxonomizr package.

### Organism(s)
primer BLAST defaults to homo sapiens, so it is important that you supply a specific organism or organisms. NCBI's taxids can be found [here](https://www.ncbi.nlm.nih.gov/taxonomy). You can specify multiple organism by passing a character vector containing each of the options, like in the example below.

## Search options

get_blast_seeds passes many parameters to NCBI's primer blast tool. You can match the parameters to the fields available in the GUI here. First, use your browser to view the page source. Search for the field you are interested in by searching for the title of the field. It should be enclosed in a tag. Inside the label tag, it says for = "<name_of_parameter>". Copy the string after for = and add it to get_blast_seeds as the name of a parameter, setting it equal to whatever you like.

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
You can find the description and suggested values for this search option. HITSIZE ='1000000' is added to the search below along with several options that increase the number of entries returened from primersearch.

```
blast_seeds_parent <- "/Users/limeybean/Dropbox/CRUX_2.0/12S_V5F1_modified_params"
metabarcode <- "12S_V5F1"
accession_taxa_path <- "/Users/limeybean/Dropbox/CRUX_2.0/accession2taxid/accessionTaxa.sql"

get_blast_seeds("TAGAACAGGCTCCTCTAG", "TTAGATACCCCACTATGC",
                 blast_seeds_parent, metabarcode, accession_taxa_path, HITSIZE ='1000000', evalue='100000',
                 word_size='6', MAX_TARGET_PER_TEMPLATE = '5',
                 NUM_TARGETS_WITH_PRIMERS ='500000', minimum_length = 50,
                 MAX_TARGET_SIZE = 200,
                 organism = c("1476529", "7776"), return_table = FALSE)

# This results in 111500 blast seed returns, note the default generated 1047.

```
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





# [rcrux_blast](https://lunagal.github.io/rcrux_blast)

rcrux_blast uses the entries generated by get_blast_seeds and the nucleotide-nucleotide matching of blastn to generate a .csv of ncbi database entries that match a sequence found in the get_blast_seeds step.

## Internal data pipeline

rcrux_blast is a wrapper function that passes data to blast_datatable. rcrux_blast handles the creation of a hidden save directory and an output directory and writes a .csv summarizing the results of blast_datatable to the output directory. Optionally, it also writes .csvs detailing entries rejected by blast_datatable.

Internally, blast_datatable repeatedly samples rows from the table of seeds, calls blastdbcmd on each accession number, and uses blastn to build a table of nucleotide matches. It samples by drawing random indices from a list of unsampled indices and examining the rows at those indices. It passes those rows to run_blastdbcmd, which extracts the accession number, forward stop, and reverse stop, then uses them as arguments for blastdbcmd. blastdbcmd outputs a fasta, which blast_datatable aggregates into a multi-fasta character vector. blast_datatable purges any entry that has more than a specified number of Ns or did not return a result, recording those indices. When it has finished building the mutli-fasta vector, it passes it to blastn, which returns every nucleotide sequence that matches a sequence in the file. run_blastn parses the blastn output into a data.frame, and blast_datatable adds that data.frame to its output. It repeats this process until it has sampled every row. Then, it uses taxonomizr to add taxonomic data to the data.frame based on the accession numbers. The final output is the aggregate of all blastn calls with the taxonomic data added.

## Example

In this example, rcrux_blast is called on the .csv generated by the get_blast_seeds example above. This example does not rely on an internet connection, but it will use a lot of memory and CPU time.

```
# These file directories need to be changed to locations on your device
blast_seeds_parent <- "D:/blast_seeds_test"
accession_taxa_path <- "D:/taxonomizr_data/accessionTaxa.sql"

# This path is indepedent of device; it only depends on get_blast_seeds
# having been run with Metabarcode_name = "12S_V5F1"
seeds_csv_path <- paste0(blast_seeds_parent, "/12S_V5F1/12S_V5F1_raw_primerTree_output.csv")
RCRUX.dev::rcrux_blast("short_test/12S_V5F1/12S_V5F1_primerTree_output_with_taxonomy.csv",
                      "blast_test_save", db_path, accession_taxa_path)
```

## Funding

The [CalCOFI](https://calcofi.org/) program provided funding support for this project.

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
