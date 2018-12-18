
# ==== preparing environment ====

# compile and load the main library
options(devtools.install.args = "--no-multiarch")
devtools::document(".")
devtools::install(".")

# load intermediate definition
source("inst/server/server-utils.R")
local_opts = list(wdir="P:/2015/65/2018-blm-replication/")

# ==== prepaparing intermediate data file ====

# prepare data fiel for dynamic and static estimations
source("inst/server/data-selection-static.r")
source("inst/server/data-selection-dynamic.r")





