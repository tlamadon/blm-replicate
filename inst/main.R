install.packages("packrat")

# ==== make sure the environment has all dependencies
packrat::status()
packrat::restore()

# load the main library OR compile and install it
if ("blmrep" %in% rownames(installed.packages())) {
  library(blmrep)
} else {
  options(devtools.install.args = "--no-multiarch")
  devtools::document(".") # generate documentation
  devtools::install(".",upgrade = "never")  # compile and install
}

# ===== setup the parameters ===== #

local_opts = list()
local_opts$dry_run = TRUE # runs all the step, but small maxiter and small number of starting values
local_opts$use_simulated_data = TRUE
dir.create("./tmp",showWarnings = FALSE)
dir.create("./tmp/data-tmp",showWarnings = FALSE)
local_opts$wdir="./tmp"

local_opts$number_of_clusters = 15    # number of cores that are available
local_opts$bootstrap_nreps    = 10  # number of replications to use for bootstrap

# ==== prepapre options for running all results =====
source("inst/server/server-utils.R")

# ===== simualte data if dry-run ===== #
if (local_opts$use_simulated_data) {
  generate_simualted_data()
}

# ==== construct intermediate data files ====

# prepare data file for dynamic and static estimations
if (!file.exists(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))) {
  source("inst/server/data-selection-static.r")
}

#if (!file.exists(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))) {
#  source("inst/server/data-selection-dynamic.r")
#}

# ===== static estimation ========
source("inst/server/estimation-static.r")
source("inst/server/fig-blm.R")

# estimate groups for static model & save descritptive statistics
server.static.d2003.clustering()         #  ~ 0.5 cpu.h
server.static.d2003.clustering.stats()   # short

# Main mixture results & bootstrap
server.static.mixture.d2003.estimate()    # ~ 12        cpu.h
server.static.mixture.estimate.boostrap() # ~ 12*nreps  cpu.h  -- this is very long!

# Generate main figure
fig.static.mixt.means()

# Generate appendix tables
table.movers.count()
table.movers.wages()

# Compute counterfactuals

# Main regression results
server.static.mini.estimate.main() # short

# Model iteration - reclassifying
server.static.mixt.estimate.dirty_iteration() # ~ 100 cpu.h

# Mixture of mixture model
server.static.mixture.mixtofmixt() # ~ 12 cpu.h
