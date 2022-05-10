# Install the tweaked version of primerTree
# The first line of code is only necessary if you have a different version of primerTree installed already
# If you already use primerTree this may interfere with scripts that use the primer_search function

remove.packages("primerTree")

devtools::install_github("LunaGal/primerTree")