install.packages("lubridate")
install.packages("XML")
install.packages("httr")
install.packages("tidyverse")
install.packages("tidyr")
install.packages("dplyr")
install.packages("ape")
install.packages("tibble")
install.packages("rlist")
install.packages("rlang")
install.packages("taxonomizr")
install.packages("data.table")
install.packages("RCurl")
install.packages("parallel")

remove.packages("primerTree")
if(!require("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
devtools::install_github("LunaGal/primerTree")

if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("ShortRead")
