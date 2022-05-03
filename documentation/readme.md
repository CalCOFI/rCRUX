# get_blast_seeds overview

get_blast_seeds() is a function that makes a request to NCBI's [primer designer tool](https://www.ncbi.nlm.nih.gov/tools/primer-blast/) then processes it into two .csv files, one with taxonomy data and one without. It also generates two .pdf files with histograms.

# setup

get_blast_seeds requires several libraries to run. Install any packages **other than primer_tree** in the libraries.R file that you don't have installed then run libraries.R.

Run the functions.R file to define get_blast_seeds and all of its ancillary functions.

get_blast_seeds uses a local sqlite database to associate taxonomic data with the data it downloads from ncbi. The accession_taxa_setup.R file contains a script to build the sqlite database. Edit line 4 to the appropriate directory then run the whole script. (It can be whatever directory works best for you, but be sure to take note of it because you will need it for get_blast_seeds.)

At this point, you are ready to use get_blast_seeds. You can run the example.R file if you like. Edit the third and fourth lines to the appropriate directories. blast_seeds_parent is the directory that the results folder will be created inside. accesssion_taxa_path is the full path to sqlite database created by accession_taxa_setup.R.

# A Note on primerTree

The current version of get_blast_seeds relies on a tweaked version of primer_search. Instead of the canonical primerTree, install the fork by running the script in tweaked_primerTree.R. If you do not have any version of primerTree installed, you can run just line 6.

# get_blast_seeds usage

get_blast_seeds(forward_primer, reverse_primer,
                file_out_dir, Metabarcode_name,
                accessionTaxa, 
                organism, mismatch = 3,
                minimum_length = 5, maximum_length = 500,
                num_permutations = 25,
                primer_specificity_database = "nt",
                hitsize='1000000', evalue='100000',
                word_size='6',
                MAX_TARGET_PER_TEMPLATE = '5000',
                NUM_TARGETS_WITH_PRIMERS ='500000', ...,
                return_table = TRUE)
                
forward_primer, reverse_primer, organism, and every argument from num_permutations to ... (inclusive) are all passed to primer_search, which in turn passes them as strings to primer blast as search parameters.

#### mismatch, minimum_length, and maximum_length are used to filter the raw data

###### mismatch

parse_primer_hits (from Primer Tree) returns a table with mismatch_forward and mismatch_reverse columns. get_blast_seeds removes each row that has a value greater than mismatch in either of those columns

###### minimum_length and maximum_length

parse_primer_hits returns a table with a product_length column. get_blast_seeds removes each row that has a value less than the minimum or greater than the maximum in the product_length column

<!--- To Do: figure out exactly what those values mean and describe them based on this
--->

#### file_out_dir and Metabarcode_name determine where files are created

file_out_dir is the parent directory to place the data in. get_blast_seeds does not attempt to create a new directory, so file_out_dir should already exist

Metabarcode_name is used to name the subdirectory and the files. if a directory named Metabarcode_name does not exist in file_out_dir, a new directory will be created. get_blast_seeds appends Metabarcode_name to the beginning of each of the four files it generates

#### return_table

By default, get_blast_seeds returns a tibble containing the processed data. If you are running it solely for its side effects, you can disable the return by setting this parameter to false.

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
