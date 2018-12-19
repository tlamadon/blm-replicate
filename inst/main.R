
# ==== preparing environment ====

# compile and load the main library
options(devtools.install.args = "--no-multiarch")
devtools::document(".")
devtools::install(".")

# load intermediate definition
local_opts = list(wdir="P:/2015/65/2018-blm-replication/")
source("inst/server/server-utils.R")

# ==== prepaparing intermediate data file ====

# prepare data fiel for dynamic and static estimations
source("inst/server/data-selection-static.r")
source("inst/server/data-selection-dynamic.r")

# ===== static estimation ========
source("inst/server/estimation-static.r")

# estimate groups for static model & save descritptive statistics
server.static.d2003.clustering()
server.static.d2003.clustering.stats()

# Main mixture results & bootstrap
server.static.mixture.d2003.estimate()

# Compute counterfactuals

# Main regression results
server.static.mini.estimate.main()
