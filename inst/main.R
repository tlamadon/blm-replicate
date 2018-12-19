
# ==== make sure the environment has all dependencies
packrat::status()

# ==== make sure the blmrep package is available ====

# load the main library OR compile and install it
if ("blmrep" %in% rownames(installed.packages())) {
  library(blmrep)
} else {
  options(devtools.install.args = "--no-multiarch")
  devtools::document(".") # generate documentation
  devtools::install(".")  # compile and install
}

# ==== prepapre options for running all results =====

local_opts = list(wdir="P:/2015/65/2018-blm-replication/")
source("inst/server/server-utils.R")

# ==== construct intermediate data files ====

# prepare data fiel for dynamic and static estimations
if (!file.exists(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))) {
  source("inst/server/data-selection-static.r")
}
#if (!file.exists(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))) {
#  source("inst/server/data-selection-dynamic.r")
#}

# ===== static estimation ========
source("inst/server/estimation-static.r")

# estimate groups for static model & save descritptive statistics
server.static.d2003.clustering()
server.static.d2003.clustering.stats()

# Main mixture results & bootstrap
server.static.mixture.d2003.estimate()
server.static.mixture.estimate.boostrap()

# Generate main figure
fig.static.mixt.means()

# Generate appendix tables
table.movers.count()
table.movers.wages()

# Compute counterfactuals

# Main regression results
server.static.mini.estimate.main()
