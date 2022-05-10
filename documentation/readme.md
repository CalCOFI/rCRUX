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