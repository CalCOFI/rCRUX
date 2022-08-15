install.packages("lubridate", repos='http://cran.us.r-project.org')
install.packages("XML", repos='http://cran.us.r-project.org')
install.packages("httr", repos='http://cran.us.r-project.org')
install.packages("tidyverse", repos='http://cran.us.r-project.org')
install.packages("tidyr", repos='http://cran.us.r-project.org')
install.packages("dplyr", repos='http://cran.us.r-project.org')
install.packages("ape", repos='http://cran.us.r-project.org')
install.packages("tibble", repos='http://cran.us.r-project.org')
install.packages("rlist", repos='http://cran.us.r-project.org')
install.packages("rlang", repos='http://cran.us.r-project.org')
install.packages("taxonomizr", repos='http://cran.us.r-project.org')
install.packages("data.table", repos='http://cran.us.r-project.org')
install.packages("RCurl", repos='http://cran.us.r-project.org')

if(!require("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
devtools::install_github("LunaGal/primerTree")

if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("ShortRead")
