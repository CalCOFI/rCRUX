# Uses the taxonomizr library to generate an sqlite database

library(taxonomizr)

# Make sure to edit this line
taxonomizr_directory <- "/data/home/galoscarleo/taxonomy"

wd <- getwd()
dir.create(taxonomizr_directory)
setwd(taxonomizr_directory)
prepareDatabase('accessionTaxa.sql')
setwd(wd)