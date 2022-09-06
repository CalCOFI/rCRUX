# rCRUX: Generate CRUX metabarcoding reference libraries in R



**Authors:** [Luna Gal](<galoscarleo@gmail.com>), [Zachary Gold](zack.gold@ucla.edu), [Emily Curd](eecurd@g.ucla.edu)<br/>
**License:**
[GPL-3](https://opensource.org/licenses/GPL-3.0)


eDNA metabarcoding is increasingly used to survey biological communities using common universal and novel genetic loci. There is a need for an easy to implement computational tool that can generate metabarcoding reference libraries for any locus, and are specific, and comprehensive. We have reimagined CRUX (Curd et al. 2019) and developed the rCRUX package to fit this need by generating taxonomy and miltifasta files for any user defined locus.  The typical data pipeline involves using get_blast_seeds() to download and wrangle results from NCBI’s primer BLAST tool, then using rcrux_blast() to search a local NCBI database for matches.


## Workflow

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

NCBI's BLAST+ suite must be locally installed. NCBI provides installation instructions for [Windows](https://www.ncbi.nlm.nih.gov/books/NBK52637/), [Linux](https://www.ncbi.nlm.nih.gov/books/NBK52640/), and [Mac OS](https://www.ncbi.nlm.nih.gov/books/NBK569861/). Additionally, it is available in the Debian archives as ncbi-blast+ and as a homebrew formula called blast.

**Blast-formatted database**

rCRUX requires a local blast-formatted database. NCBI provides a tool for downloading databases as part of the blast+ package. A brief help page can be found [here](https://www.ncbi.nlm.nih.gov/books/NBK569850/).

**Taxonomizr**

rCRUX uses the taxonomizr package for taxonomic assignment based on NCBI [Taxonomy id's \(taxids\)](https://www.ncbi.nlm.nih.gov/). taxonomizr requires a local sqlite database and provides prepareDatabase to automatically build the local database. Many RCRUX functions require a path to a taxonomizr database, so you should run prepareDatabase before running RCRUX functions. You can read about taxonomizr [here](https://cran.r-project.org/web/packages/taxonomizr/vignettes/usage.html).


# Example pipeline

The following example shows a simple RCRUX pipeline from start to finish. Note that this example will require internet access and considerable database storage, run time (mainly for blastn), and system resources to execute.


```
seeds_path <- "/my/rCRUX_output_directory"
accession_path <- "/my/accessionTaxa.sql"
blastdb_path <- "/my/local/blast_database/nt"

get_blast_seeds("TAGAACAGGCTCCTCTAG", "TTAGATACCCCACTATGC",
                 blast_seeds_parent, "12S_V5F1", accession_path,
                 organism = c("7776", "7777"), return_table = FALSE)


# A .csv is automatically created at this path based on the arguments passed to get_blast_seeds

csv_path <- "/my/directory/12S_V5F1/12S_V5F1_primerTree_output_with_taxonomy.csv"

rcrux_blast(csv_path, "blast_test_save", db_path, accession_taxa_path)

```




# [get_blast_seeds](https://lunagal.github.io/get_blast_seeds)

This script takes a forward and reverse primer sequence and using get_blast_seeds() generates a csv with data returned from [NCBI's primer blast](https://www.ncbi.nlm.nih.gov/tools/primer-blast/) about full length barcode sequence containing primer matches.

get_blast_seeds uses modified versions of functions from the primerTree package to submit queries to NCBI's primer BLAST tool, then aggregates results into a single data.frame. primer_search expands degenerate primers into each possible non-degenerate primer and submits a query for each. get_blast_seeds further multiplies the number of queries by allowing the user to query the primers for each organism in a vector. get_blast_seeds collects all these results from primer_search, filters them based on product length, and adds taxonomic data using the taxonomizr package.


# get_blast_seeds overview

get_blast_seeds() is a function that makes a request to NCBI's [primer designer tool](https://www.ncbi.nlm.nih.gov/tools/primer-blast/) then processes it into two .csv files, one with taxonomy data and one without. It also generates two .pdf files with histograms.

# setup

get_blast_seeds requires several libraries to run. Install each of the following packages that you do not yet have installed:

```
install.packages(lubridate)
install.packages(XML)
install.packages(httr)
install.packages(tidyverse)
install.packages(tidyr)
install.packages(dplyr)
install.packages(ape)
install.packages(tibble)
install.packages(rlist)
install.packages(rlang)
install.packages(taxonomizr)
install.packages(data.table)
install.packages(RCurl)
install.packages(parallel)
```
Additionally, run the following code block to install a modified version of primerTree. See the seciont below for more information.

```
# If you already use primerTree this may interfere with scripts that use the primer_search function

remove.packages("primerTree")

devtools::install_github("LunaGal/primerTree")
```
Run the functions.R file to define get_blast_seeds and all of its ancillary functions.

get_blast_seeds uses a local sqlite database to associate taxonomic data with the data it downloads from ncbi. The accession_taxa_setup.R file contains a script to build the sqlite database. Edit line 4 to the appropriate directory then run the whole script. (It can be whatever directory works best for you, but be sure to take note of it because you will need it for get_blast_seeds.)

```
# Uses the taxonomizr library to generate an sqlite database

library(taxonomizr)

# Make sure to edit this line
taxonomizr_directory <- "D:/test"

wd <- getwd()
dir.create(taxonomizr_directory)
setwd(taxonomizr_directory)
prepareDatabase('accessionTaxa.sql')
setwd(wd)
```

At this point, you are ready to use get_blast_seeds. You can run the fragment below, taken from the example.R file if you like. Edit the third and fourth lines to the appropriate directories. blast_seeds_parent is the directory that the results folder will be created inside. accesssion_taxa_path is the full path to sqlite database created by accession_taxa_setup.R.

```
# These file directories need to be changed to locations on your device

blast_seeds_parent <- "D:/blast_seeds"
accession_taxa_path <- "D:/test/accessionTaxa.sql"

get_blast_seeds("TAGAACAGGCTCCTCTAG", "TTAGATACCCCACTATGC",
                 blast_seeds_parent, "12S_V5F1", accession_taxa_path,
                 organism = c("7776", "7777"), return_table = FALSE)
```

# A Note on primerTree

The current version of get_blast_seeds relies on a tweaked version of primer_search. Instead of the canonical primerTree, install the fork by running the script in tweaked_primerTree.R. If you do not have any version of primerTree installed, you can run just line 6. Note that because this slightly changes the behavior of primer_search, it may interfere with any scripts that use the canonical version of primerTree.

# get_blast_seeds usage

get_blast_seeds(forward_primer,
                reverse_primer,
                file_out_dir,
                Metabarcode_name,
                accessionTaxa,
                organism, mismatch = 3,
                minimum_length = 5,
                maximum_length = 500,
                num_permutations = 25,
                primer_specificity_database = "nt",
                hitsize='1000000',
                evalue='100000',
                word_size='6',
                MAX_TARGET_PER_TEMPLATE = '5000',
                NUM_TARGETS_WITH_PRIMERS ='500000',
                ...,
                return_table = TRUE)

#### Arguments

| Argument | Description  |
| -- | -- |
| ``forward_primer`` | passed to primer_search, which turns it into a list of each primer it could be based on its degenerate primers, then passes each one in turn to NCBI |
| ``reverse_primer`` | passed to primer_search, which turns it into a list of each primer it could be based on its degenerate primers, then passes each one in turn to NCBI |
| ``file_out_dir`` | the parent directory to place the data in. get_blast_seeds does not attempt to create a new directory, so file_out_dir should already exist |
| ``Metabarcode_name`` | used to name the subdirectory and the files. If a directory named Metabarcode_name does not exist in file_out_dir, a new directory will be created. get_blast_seeds appends Metabarcode_name to the beginning of each of the four files it generates. |
| ``accessionTaxa`` | the path to sql created by taxonomizr |
| ``organism`` | a vector of character vectors. Each character vector is passed in turn to primer_search, which passes them to NCBI. get_blast_seeds aggregates all of the results into a single file. |
| ``minimum_length`` | parse_primer_hits returns a table with a product_length column. get_blast_seeds removes each row that has a value less than minimum_length in the product_length column |
| ``maximum_length`` | parse_primer_hits returns a table with a product_length column. get_blast_seeds removes each row that has a value greater than maximum_length in the product_length column |
| ``num_permutations`` | passed to primer_search, which passes it to NCBI |
| ``primer_specificity_database`` | passed to primer_search, which passes it to NCBI |
| ``hitsize`` | passed to primer_search, which passes it to NCBI |
| ``evalue`` | passed to primer_search, which passes it to NCBI |
| ``word_size`` | passed to primer_search, which passes it to NCBI |
| ``MAX_TARGET_PER_TEMPLATE`` | passed to primer_search, which passes it to NCBI |
| ``NUM_TARGETS_WITH_PRIMERS`` | passed to primer_search, which passes it to NCBI |
| ``...`` | additional arguments passed to primer_search, which passes it to NCBI |
| ``return_table`` | determines whether get_blast_seeds returns a tibble containing the same information as the csv it creates |

#### Example

```
# These file directories need to be changed to locations on your device

blast_seeds_parent <- "D:/blast_seeds"
accession_taxa_path <- "D:/test/accessionTaxa.sql"

get_blast_seeds("TAGAACAGGCTCCTCTAG", "TTAGATACCCCACTATGC",
                 blast_seeds_parent, "12S_V5F1", accession_taxa_path,
                 organism = c("7776", "7777"), return_table = FALSE)
```

An example output obtained using these parameters at 6:05 EST 9 May 2022 can be found in examples.

# NCBI primer blast parameters

As discussed above, get_blast_seeds passes many parameters to NCBI's primer blast tool. You can match the parameters to the fields available in the GUI [here](https://www.ncbi.nlm.nih.gov/tools/primer-blast/). First, use your browser to view the page source. Search for the field you are interested in by searching for the title of the field. It should be enclosed in a <label> tag. Inside the label tag, it says for = "[name_of_parameter]". Copy the string after for = and add it to get_blast_seeds as the name of a parameter, setting it equal to whatever you like.

Example: I want to set "Exon junction span" to 10. I open the source of the primer designing tool and look for that string. I find the following:

```<label class="m" for="PRIMER_ON_SPLICE_SITE">Exon junction span</label>```

I copy PRIMER_ON_SPLICE_SITE and add it to get_blast_seeds:

```
get_blast_seeds("TAGAACAGGCTCCTCTAG", "TTAGATACCCCACTATGC",
                blast_seeds_parent, "12S_V5F1", accession_taxa_path,
                organism = c("7776"), MAX_TARGET_PER_TEMPLATE = 10,
                PRIMER_ON_SPLICE_SITE = "10"
                return_table = FALSE)
```

# Selected common errors

#### _ is not a valid url. It will be ignored.

This means that NCBI returned a response that parse_primer_hits couldn't parse. Usually, this is because the query was too large or had an argument that NCBI didn't know how to handle, such as an organism that doesn't exist. Pasting the URL into a browser can sometimes help determine which. It may have an error message that reads, in part

```
Invalid organism or taxonomy id input: dogpf.
```

If your query was very large, rerunning it later and hoping NCBI is less trafficked at that time may help.

#### Error in accessionToTaxa(input$accession, accessionTaxa_path) : _ does not exist.

This means either that you did not enter the correct path for the accessionTaxa argument. You may not have edited the line ``accession_taxa_path <- "xyz"`` (and if you didn't make sure you also edit the previous line) or you may not have successfully created the accessionTaxa sql.




## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Epskamp2018ggm" class="csl-entry">

Epskamp, Sacha, Lourens J. Waldorp, Rene Mottus, and Denny Borsboom.
2018. “<span class="nocase">The Gaussian Graphical Model in
Cross-Sectional and Time-Series Data</span>.” *Multivariate Behavioral
Research* 53 (4): 453–80.
<https://doi.org/10.1080/00273171.2018.1454823>.

</div>

<div id="ref-williams2019nonregularized" class="csl-entry">

Williams, Donald R., Mijke Rhemtulla, Anna C Wysocki, and Philippe Rast.
2019. “On Nonregularized Estimation of Psychological Networks.”
*Multivariate Behavioral Research* 54 (5): 719–50.
<https://doi.org/10.1080/00273171.2019.1575716>.

</div>

</div>
