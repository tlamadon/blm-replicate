
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

# estimate groups for static model
server.static.d2003.clustering()
# save group stats
server.static.d2003.clustering.stats()


# Mini model
server.static.mini.estimate.main()
