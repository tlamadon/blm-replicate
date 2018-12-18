
# ==== preparing environment ====

# compile and load the main library
devtools::document(".")
devtools::install(".")

# load intermediate definition
source("inst/server/server-utils.R")


# ==== prepaparing intermediate data file ====

# prepare data fiel for dynamic and static estimations
source("inst/server/data-selection-static.r")
source("inst/server/data-selection-dynamic.r")





