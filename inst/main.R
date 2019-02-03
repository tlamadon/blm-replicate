# ===== setup the parameters ===== #
local_opts = list()
dir.create("./tmp",showWarnings = FALSE)
dir.create("./tmp/data-tmp",showWarnings = FALSE)
local_opts$wdir="./tmp"

# ---------   dry run configuration ---------
local_opts$use_simulated_data = TRUE
local_opts$cpu_count                 = 2    # number of cores that are available
local_opts$bootstrap_nreps           = 10  # number of replications to use for bootstrap
local_opts$estimation.mixture = list(maxiter = 50,est_rep=4,est_nbest=2)

# ==== prepapre options for running all results =====
source("inst/server/server-utils.R")
source("inst/server/estimation-static.r")
source("inst/server/fig-blm.R")
generate_simualted_data(force=TRUE);

# ==== construct intermediate data files ====

# prepare data file for dynamic and static estimations
if (!file.exists(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))) {
  source("inst/server/data-selection-static.r")
}

#if (!file.exists(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))) {
#  source("inst/server/data-selection-dynamic.r")
#}

# ===== static estimation ========

# estimate groups for static model & save descritptive statistics
server.static.d2003.clustering()         #  ~ 0.5 cpu.h
server.static.d2003.clustering.stats()   # short

# Main mixture results & bootstrap
server.static.mixture.d2003.estimate()    # ~ 12        cpu.h
server.static.mixture.estimate.boostrap() # ~ 12*nreps  cpu.h  -- this is very long!

# Generate main figure
fig.static.mixt.means()

# Generate appendix tables & figures
tab.static.movers.count()
fig.static.movers.wages()
fig.dynamic.mixt.connectedness() # connectedness picture

# Compute counterfactuals
server.static.analysis.meaneffects()

# Create the main table
tab.static.mixt.vdec()

# Main regression results
server.static.mini.estimate.main()      # short
server.static.mini.estimate.bootstrap() # ~ 2 cpu.h

# Model iteration - reclassifying
server.static.mixt.estimate.model_iteration() # ~ 100 cpu.h

# Robustness
server.static.mixture.mixtofmixt()            #  Mixture of mixture model ~ 12 cpu.h
server.static.mixture.estimate.robust.fsize() # less than 50, larger than 50
server.static.mixture.estimate.robust.nf()    # varying number of firm types
server.static.mixture.estimate.robust.nk()    # varying number of worker types

# ===== dynamic estimation ========

table.statedependence()
table.endogeneousMobility()


# ====== probabilistic estimation ========


# ====== shimer-smith simulation and estimation =======





# ======= preparation =======
# if (!"packrat" %in% rownames(installed.packages())) {
#   install.packages("packrat")
# }
#
# # ==== make sure the environment has all dependencies
# packrat::set_opts(ignored.packages = c("blmrep"))
# packrat::status()
# packrat::restore()
#
# # load the main library OR compile and install it
# if ("blmrep" %in% rownames(installed.packages())) {
#   library(blmrep)
# } else {
#   options(devtools.install.args = "--no-multiarch")
#   #devtools::document(".") # generate documentation
#   devtools::install(".",upgrade = "never")  # compile and install
# }
