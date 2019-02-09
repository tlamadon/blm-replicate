require(blmrep)

# ===== setup the parameters ===== #
local_opts = list()
dir.create("./tmp",showWarnings = FALSE)
dir.create("./tmp/data-tmp",showWarnings = FALSE)
local_opts$wdir="./tmp"

local_opts$use_simulated_data       = FALSE
local_opts$cpu_count                = detectCores()    # number of cores that are available
local_opts$bootstrap_nreps          = 200  # number of replications to use for bootstrap
local_opts$estimation.mixture       = list(maxiter = 2000,est_rep=50,est_nbest=10)
local_opts$estimation.probabilistic = list(maxiter = 200,gibbs_nfirm=100)
local_opts$trace$nm_list            = c(1,2,3,4,5,10,15,20,50,100)

# ==== prepapre options for running all results =====
source("inst/server/server-utils.R")
source("inst/server/estimation-static.r")
source("inst/server/estimation-dynamic.r")
source("inst/server/fig-blm.R")
generate_simulated_data();

# ==== construct intermediate data files ====

# prepare data file for dynamic and static estimations
if (!file.exists(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))) {
  source("inst/server/data-selection-static.r")
}
if (!file.exists(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))) {
  source("inst/server/data-selection-dynamic.r")
}

# ===== static estimation ========

# estimate groups for static model & save descritptive statistics
server.static.d2003.clustering()         #  ~ 0.5 cpu.h
server.static.d2003.clustering.stats()   # short
table.static.clusters()

# Main mixture results & bootstrap
server.static.mixture.d2003.estimate()    # ~ 12        cpu.h
server.static.mixture.estimate.boostrap() # ~ 12*nreps  cpu.h  -- this is very long!

# Generate main figure
fig.static.mixt.means()

# Generate appendix tables & figures
tab.static.movers.count()
fig.static.movers.wages()
fig.static.mixt.connectedness() # connectedness picture

# Compute counterfactuals
server.static.analysis.meaneffects()

# Create the main table
tab.static.mixt.vdec()

# Model iteration - reclassifying
server.static.mixt.estimate.model_iteration() # ~ 100 cpu.h

# Robustness
server.static.mini.estimate.main()            # linear and interacted regressions
server.static.mixture.mixtofmixt()            # Mixture of mixture model ~ 12 cpu.h
server.static.mixture.estimate.robust.fsize() # less than 50, larger than 50
server.static.mixture.estimate.robust.nf()    # varying number of firm types
server.static.mixture.estimate.robust.nk()    # varying number of worker types
server.static.mixt.estimate.robustess.residuals()
server.static.estimate.clustersplits("rk-prank")
server.static.estimate.clustersplits("rk-va")

tab.satic.robust()

# ===== dynamic estimation ========

server.dynamic.d2003.computeclusters()
server.dynamic.d2003.clustering.stats()
server.dynamic.rho.analysis()
fig.dynamic.rho_slices()

server.dynamic.mixture.d2003.estimate()
server.dynamic.mixture.d2003.boostrap()

fig.dynamic.mixt.means()
fig.dynamic.mixt.connectedness()
tab.dynamic.parameters()

# create counterfactuals
server.dynamic.analysis.meaneffects()
server.dynamic.analysis.endogenousmobility()
server.dynamic.analysis.stateDependence()

tab.dynamic.mixt.vdec()
tab.dynamic.statedependence()
tab.dynamic.endogeneousMobility()

# robustness
server.dynamic.mini.estimate()
server.dynamic.mixture.estimate.robust.nf()
server.dynamic.mixture.estimate.robust.nk()
server.dynamic.mixture.estimate.robust.different_rho()
server.dynamic.mixture.d2003.estimate.model_iteration()

tab.dynamic.robust()

# ===== Andrews, Kline, BLM comparaison ========
# server.fe.trace()
# figure.hybrid()


# ====== probabilistic estimation ========
server.static.proba.results()
fig.static.proba.gibbs()

# ====== shimer-smith simulation and estimation =======
server.shimersmith.results()
fig.shimersmith.model()
fig.shimersmith.CK_event_study()
fig.shimersmith.wages()

