# ========= UTILS =========

#' loads static data
server.static.data <- function(remove_movers_from_sdata=T) {
  load(sprintf("%s/data-tmp/data-static.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]
  if (remove_movers_from_sdata) {
    sdata=sdata[move==FALSE]
  }
  sim = list(jdata=jdata,sdata=sdata)
  rm(sdata,jdata)
  return(sim)
}

# ========= CLUSTERING =========

#' main estimation for the mixture model, this goes in
#' the paper main results
server.static.d2003.clustering <- function() {
  sim = server.static.data()

  # we start with clustering
  set.seed(65542134)
  ms    = grouping.getMeasures(sim,"ecdf",Nw=20,y_var = "y1")
  grps  = grouping.classify.once(ms,k = 10,nstart = 10000,iter.max = 200,step=100)
  sim   = grouping.append(sim,grps$best_cluster)
  res.save("m2-mixt-d2003-groups",grps)
}

server.static.d2003.clustering.stats <- function() {

  # computing some statistics
  sim = server.static.data(remove_movers_from_sdata=FALSE)

  # get clsuters
  grps  = res.load("m2-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  mstats = sim$jdata[,list(m1=mean(y1),sd1=sd(y1),
                           m2=mean(y2),sd2=sd(y2),
                           v12 = cov(y1,y2),.N),list(j1,j2)]
  cstats = sim$sdata[,list(m1=mean(y1),sd1=sd(y1),
                           m2=mean(y2),sd2=sd(y2),
                           v12 = cov(y1,y2),.N),list(j1)]

  sdata = sim$sdata
  setkey(sdata,j1)
  sdata[,age:=2002 - birthyear]
  sdata[,asize := .N,f1]
  sdata[,fmw := mean(y1),f1]

  gstats = rBind(
    sdata[, data.frame(get.stats.clusters(.SD,ydep="y1")), j1],
    data.frame(list(get.stats.clusters(sdata,ydep="y1"),j1=0)))

  res.save("m2-stats-d2003",list(mstats=mstats,cstats=cstats,gstats = gstats))
}

# ========== MINI-MODEL ===========

server.static.mini.d2003.estimate <- function() {

  sim = server.static.data()
  # get clsuters
  grps  = res.load("m2-mixt-y2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  # --- MINI ESTIMATION ---- #
  set.seed(5422349)
  rs         = rkiv0.start("m2-mini-d2003-profandlinear")
  rs$info    = "static mini model on 2003 data, using profiling and with B=1"

  model_mini = m2.mini.estimate(sim$jdata,sim$sdata,method="prof")
  model_minilin      = m2.mini.estimate(sim$jdata,sim$sdata,method="linear")
  rkiv0.put(rs,list(mini=model_mini,mini_lin=model_minilin))
}


server.static.mini.estimate.main <- function(){
  sim = server.static.data()
  grps  = res.load("m2-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  mini_model = m2.mini.estimate(sim$jdata,sim$sdata,norm = 1,method="prof")
  res.save("m2-mini-prof",mini_model)

  mini_model = m2.mini.estimate(sim$jdata,sim$sdata,norm = 1,method="linear")
  res.save("m2-mini-linear",mini_model)
}



# ==== MIXTURE MODEL =====

server.static.mixture.d2003.estimate <- function() {

  sim = server.static.data()
  grps  = res.load("m2-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  # --- MAIN ESTIMATION ---- #
  set.seed(87954352)

  # --- use cluster ---#
  cl = makeCluster(local_opts$cpu_count)
  clusterEvalQ(cl,require(blmrep))

  # --- MAIN ESTIMATION -  STATIONARY INTERACTIONS ---- #
  ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=local_opts$estimation.mixture$est_rep,
                         est_nbest=local_opts$estimation.mixture$est_nbest,
                         sdata_subsample=0.1,
                         maxiter = local_opts$estimation.mixture$maxiter)

  res_mixt = m2.mixt.estimate.all(sim,nk=6,ctrl,cl)
  res.save("m2-mixt-d2003-main-fixb",res_mixt)

  # --- MAIN ESTIMATION - NON STATIONARY ---- #
  ctrl$fixb=FALSE
  res_mixt = m2.mixt.estimate.all(sim,nk=6,ctrl,cl)
  res.save("m2-mixt-d2003-main-ns",res_mixt)

  stopCluster(cl)
}


server.static.mixture.d2003.fit <- function() {

  load(sprintf("%s/data-tmp/data-static.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]
  sim = list(jdata=jdata,sdata=sdata)
  rm(sdata,jdata)

  # get clsuters
  grps  = res.load("m2-mixt-y2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  # get main estimate
  res_main = res.load("m2-mixt-y2003-main-fixb")

  # impute data on movers
  jdata.sim = m2.mixt.impute.movers(sim$jdata,res_main$model)[,list(y1,y1_imp,j1,j2,y2,y2_imp)]

  vmin = quantile(jdata.sim$y1,0.01)
  vmax = quantile(jdata.sim$y1,0.99)

  jdata.sim[, j1b := ceiling(j1/2)]
  jdata.sim[, j2b := ceiling(j2/2)]

  dd_m = jdata.sim[, {
    d1 = density(y1,bw="sj",from=vmin,to=vmax,n=30)
    d2 = density(y1_imp,bw="sj",from=vmin,to=vmax,n=30)
    data.frame(x=d1$x,y1=d1$y,y2=d2$y)
  }, list(j1b=sprintf("%i-%i",2*j1b-1,2*j1b),j2b=sprintf("%i-%i",2*j2b-1,2*j2b)) ]

  rr_m = m2.movers.checkfit(jdata.sim)

  sdata.sim = m2.mixt.impute.stayers(sim$sdata,res_main$model)
  rr_s      = m2.stayers.checkfit(sdata.sim)

  vmin = quantile(sdata.sim$y1,0.001)
  vmax = quantile(sdata.sim$y1,0.999)

  dd_s = sdata.sim[, {
    d1 = density(y1,from=vmin,to=vmax,n=30)
    d2 = density(y1_imp,from=vmin,to=vmax,n=30)
    data.frame(x=d1$x,y1=d1$y,y2=d2$y)
  }, list(j1) ]

  rs    = rkiv0.start("m2-mixt-d2003-main-fit",info="fit of the model")
  rkiv0.put(rs,list(dd_s=dd_s,dd_m=dd_m,rr_s=rr_s,rr_m=rr_m))
}

#' this estimates the model with different values of K to check
#' robustness. We recluster every time
server.static.mixture.estimate.robust.nf <- function() {

  sim = server.static.data()
  grps  = res.load("m2-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  set.seed(87954352)
  ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=local_opts$estimation.mixture$est_rep,
                         est_nbest=local_opts$estimation.mixture$est_nbest,
                         sdata_subsample=0.1,
                         maxiter = local_opts$estimation.mixture$maxiter)

  cl = makeCluster(local_opts$cpu_count)
  clusterEvalQ(cl,require(blmrep))

  rr_mixt = list()
  for (nf_size in c(3,5,20)) {
    tryCatch({
      ms    = grouping.getMeasures(sim,"ecdf",Nw=20,y_var = "y1")
      grps  = grouping.classify.once(ms,k = nf_size,nstart = 500,iter.max = 200)
      sim   = grouping.append(sim,grps$best_cluster)
      res_mixt = m2.mixt.estimate.all(sim,nk=6,ctrl,cl=cl)
      rr_mixt[[paste("nf",nf_size,sep="-")]] = res_mixt
    })}
  stopCluster(cl)

  res.save("m2-mixt-d2003-change-nf",rr_mixt)
}

#' this estimates the model with different values of K to check
#' robustness. We recluster every time
server.static.mixture.estimate.robust.nk <- function() {

  sim = server.static.data()
  grps  = res.load("m2-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=local_opts$estimation.mixture$est_rep,
                         est_nbest=local_opts$estimation.mixture$est_nbest,
                         sdata_subsample=0.1,
                         maxiter = local_opts$estimation.mixture$maxiter)

  cl = makeCluster(local_opts$cpu_count)
  clusterEvalQ(cl,require(blmrep))

  rr_mixt = list()
  for (nk_size in c(3,5,9)) {
    tryCatch({
      res_mixt = m2.mixt.estimate.all(sim,nk=nk_size,ctrl=ctrl,cl=cl)
      res_mixt$second_stage_reps_all=NA
      rr_mixt[[paste("nk",nk_size,sep="-")]] = res_mixt
    })}
  stopCluster(cl)

  res.save("m2-mixt-d2003-change-nk",rr_mixt)
}


server.static.mixture.estimate.robust.fsize <-function() {
  res_mixt_lt_50 = server.static.mixture.estimate.robust.fsize.int(0,50)
  res_mixt_gt_50 = server.static.mixture.estimate.robust.fsize.int(50,Inf)
  res.save("m2-mixt-d2003-firmsize",list(mixt_leq_50=res_mixt_lt_50,mixt_gt_50=res_mixt_gt_50))
}

#' we estimate for small firms and for large firms
server.static.mixture.estimate.robust.fsize.int <- function(min_size,max_size) {

  sim = server.static.data()
  grps  = res.load("m2-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  # select firms larger or smaller than 50
  fids = unique(sim$sdata[,.N,f1][ (N>min_size) & (N<max_size) , f1])

  sim$sdata = sim$sdata[f1 %in% fids]
  sim$jdata = sim$jdata[f1 %in% fids][f2 %in% fids]

  set.seed(87954352)
  ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=local_opts$estimation.mixture$est_rep,
                         est_nbest=local_opts$estimation.mixture$est_nbest,
                         sdata_subsample=0.1,
                         maxiter = local_opts$estimation.mixture$maxiter)

  cl = makeCluster(local_opts$cpu_count)
  clusterEvalQ(cl,require(blmrep))
  ms    = grouping.getMeasures(sim,"ecdf",Nw=20,y_var = "y1")
  grps  = grouping.classify.once(ms,k = 10,nstart = 1000,iter.max = 200)
  sim   = grouping.append(sim,grps$best_cluster)
  res_mixt = m2.mixt.estimate.all(sim,nk=6,ctrl,cl=cl)
  stopCluster(cl)

  return(res_mixt)
}



#' Boostrapping the main results. We use the value for the bias-corrected
#' estimates.
server.static.mixture.estimate.boostrap <- function(){

  sim = server.static.data()
  grps  = res.load("m2-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  # get main estimate
  res_main  = res.load("m2-mixt-d2003-main-fixb")
  model_true = res_main$model

  ctrl      = em.control(nplot=1000,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=local_opts$estimation.mixture$est_rep,
                         est_nbest=local_opts$estimation.mixture$est_nbest,
                         sdata_subsample=0.1,
                         maxiter = local_opts$estimation.mixture$maxiter)

  cl = makeCluster(local_opts$cpu_count)
  clusterEvalQ(cl,require(blmrep))

  # now we bootstrap with clustering
  rrr=data.frame()
  rr_mixt = list()
  for (i in 1:local_opts$bootstrap_nreps) {
    tryCatch({
    sdata.sim = m2.mixt.impute.stayers(sim$sdata,model_true)
    jdata.sim = m2.mixt.impute.movers(sim$jdata,model_true)
    sdata.sim[,y1:=y1_imp][,y2:=y2_imp][,jt:=j1][,y1_bu:=y1]
    jdata.sim[,y1:=y1_imp][,y2:=y2_imp]
    sim.sp = list(jdata=jdata.sim,sdata=sdata.sim)

    minib    = m2.mini.estimate(jdata.sim,sdata.sim,method = "prof")
    minilinb = m2.mini.estimate(jdata.sim,sdata.sim,method = "linear")

    # --- recluster ---- #
    ms    = grouping.getMeasures(sim.sp,"ecdf",Nw=20,y_var = "y1")
    grps  = grouping.classify.once(ms,k = 10,nstart = 1000,iter.max = 200)
    sim.sp   = grouping.append(sim.sp,grps$best_cluster)

    mini    = m2.mini.estimate(sim.sp$jdata,sim$sdata,method = "prof")
    minilin = m2.mini.estimate(sim.sp$jdata,sim$sdata,method = "linear")

    # --- estimate model --- #
    model_bis            = m2.mixt.estimate.all(sim.sp,nk=6,ctrl,cl)
    model_bis$clus_cor  = sim.sp$sdata[,cor(jt,j1)]

    rr_mixt[[paste(i)]] = model_bis

    rr2 = rbind(
      data.frame(minib$vdec$stats ,where="mini-truecluster"),
      data.frame(minilinb$vdec$stats ,where="linear-truecluster"),
      data.frame(mini$vdec$stats ,where="mini"),
      data.frame(minilin$vdec$stats ,where="linear"),
      data.frame(model_bis$vdec$stats ,where="mixt"))
    rr2$rep=i
    rrr = rbind(rrr,rr2)

    flog.info("done with rep %i",i)

    print(rrr)
  }, error = function(e) {catf("error in boot strap rep %i!\n",i);print(e)})
  }

  stopCluster(cl)

  rr_mixt2 = lapply(rr_mixt,function(r) { r$second_stage_reps_all=NULL;r})
  res.save("m2-mixt-d2003-bootstrap",list(vdecs=rrr,mixt_all=rr_mixt2))
}

server.static.estimate.clustersplits <- function(measure="rk-prank") {

  sim = server.static.data()
  grps  = res.load("m2-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  frank <- sim$sdata[,mean(y1),f1][,list(fq=rank(V1)/.N,f1)]
  setkey(frank,f1)
  setkey(sim$jdata,f1)
  sim$jdata[,fq1 := frank[sim$jdata,fq]]
  setkey(sim$jdata,f2)
  sim$jdata[,fq2 := frank[sim$jdata,fq]]

  # create the measure
  # for each firm we count number of movers coming/going going to each firm-quartile
  if (measure=="rk-prank") {
    #rr = sim$jdata[,list(r=median(fq1)),list(f=f2)]
    rr = sim$jdata[,list(r=.N),list(f=f1)]
    rr2 = sim$sdata[,list(r2=.N),list(f=f1)]
    setkey(rr,f)
    setkey(rr2,f)
    rr = rr2[rr]
    rr[,r:=r/r2]
    rr[,r2:=NULL]
  } else if (measure=="rk-va") {
    rr = unique(sim$sdata[,list(r=va1),list(f=f1)])
  }

  # attach this to the data with cluster
  grps  = res.load("m2-mixt-d2003-groups")
  clus  = data.table(f = names(grps$best_cluster), j = grps$best_cluster)
  setkey(rr,f)
  setkey(clus,f)
  rr = clus[rr]

  # split in 2 within
  rr[,g:=r<median(r),j]
  rr[,jn := j + 10*g]
  clus = rr[,jn]
  names(clus) = rr[,f]

  sim = server.static.data()
  sim = grouping.append(sim,clus,drop=T)

  # we then estimate the model
  ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=local_opts$estimation.mixture$est_rep,
                         est_nbest=local_opts$estimation.mixture$est_nbest,
                         sdata_subsample=0.1,
                         maxiter = local_opts$estimation.mixture$maxiter)

  cl = makeCluster(local_opts$cpu_count)
  clusterEvalQ(cl,require(blmrep))
  res_mixt = m2.mixt.estimate.all(list(sdata=sim$sdata[move==FALSE],jdata=sim$jdata),nk=6,ctrl,cl)
  stopCluster(cl)

  res.save(sprintf("m2-mixt-d2003-clus_split-%s",measure),res_mixt)
}

#' Using the model to reclassify
server.static.mixt.estimate.model_iteration <- function() {

  # --- use cluster ---#
  cl = makeCluster(local_opts$cpu_count)
  clusterEvalQ(cl,require(blmrep))

  # --- MAIN ESTIMATION -  STATIONARY INTERACTIONS ---- #
  ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=local_opts$estimation.mixture$est_rep,
                         est_nbest=local_opts$estimation.mixture$est_nbest,
                         sdata_subsample=0.1,
                         maxiter = local_opts$estimation.mixture$maxiter)

  # ===== POACHING RANK ========= #
  # sim = server.static.data()
  # load("L:\\Tibo\\qtrdata\\tmp-2003-prank.dat")
  # fmean = data_prank[from_j2j>0][from_u>0][,list(from_u/(from_j2j + from_u),N=(from_j2j + from_u)),list(f1 = fid)]
  # clusters   = Ckmeans.1d.dp(fmean$V1, 10,fmean$N)
  # grp_akm    = clusters$cluster;names(grp_akm) = fmean$f1
  # sim        = grouping.append(sim,grp_akm,drop=T)
  #
  # res =  m2.mixt.estimate.reclassify(sim,10,ctrl,cl)
  # res.save("m2-mixt-d2003-main-reclassify-startwithprank",res)

  # ===== RATIO OF MOVERS ========= #
  sim = server.static.data(remove_movers_from_sdata = F)
  fmean = sim$sdata[,list(V1=sum(move==TRUE)/.N,.N),f1]
  sim$sdata = sim$sdata[move==FALSE]
  clusters   = Ckmeans.1d.dp(fmean$V1, 10,fmean$N)
  grp_akm    = clusters$cluster;names(grp_akm) = fmean$f1
  sim        = grouping.append(sim,grp_akm,drop=T)

  res =  m2.mixt.estimate.reclassify(sim,10,ctrl=ctrl,cl=cl)
  res.save("m2-mixt-d2003-main-reclassify-startwithratiomover",res)

  # ===== Using residuals ========= #
  # see section on residuals

  # ===== START WITH MAIN RESULTS ========= #
  sim = server.static.data()
  grps  = res.load("m2-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)
  res =  m2.mixt.estimate.reclassify(sim,10,ctrl=ctrl,cl=cl)
  res.save("m2-mixt-d2003-main-reclassify-mainresults",res)

  # ===== START WITH VALUE ADDED ========= #
  sim = server.static.data()
  fmean = sim$sdata[,log(va1[1]/size1[1]),f1][is.finite(V1)]
  clusters   = Ckmeans.1d.dp(fmean$V1, 10)
  grp_akm        = clusters$cluster
  names(grp_akm) = fmean$f1
  sim        = grouping.append(sim,grp_akm,drop=T)

  res =  m2.mixt.estimate.reclassify(sim,10,ctrl=ctrl,cl=cl)
  res.save("m2-mixt-d2003-main-reclassify-startwithva",res)

  # ===== START WITH MEANS ========= #
  sim   = server.static.data()
  fmean = sim$sdata[,list(mean(y1),.N),f1]
  clusters   = Ckmeans.1d.dp(fmean$V1, 10,fmean$N)
  grp_akm    = clusters$cluster; names(grp_akm) = fmean$f1
  sim        = grouping.append(sim,grp_akm,drop=T)

  res =  m2.mixt.estimate.reclassify(sim,10,ctrl=ctrl,cl=cl)
  res.save("m2-mixt-d2003-main-reclassify-startwithmean",res)

  # ===== START WITH AKM ========= #
  sim   = server.static.data()
  firm.fe = m2.fe.firms(sim$jdata)                 # use AKM grouping
  clusters   = Ckmeans.1d.dp(firm.fe$fe$psi, 10)   # classify using AKM FE
  grp_akm        = clusters$cluster
  names(grp_akm) = firm.fe$fe$f1
  sim        = grouping.append(sim,grp_akm,drop=T)

  res =  m2.mixt.estimate.reclassify(sim,10,ctrl=ctrl,cl=cl)
  res.save("m2-mixt-d2003-main-reclassify-startwithakm",res)

  stopCluster(cl)

  # =======  combining results ======= #
  rr_all = data.frame()
  ll = list(static="m2-mixt-d2003-main-reclassify-mainresults",
            akm   ="m2-mixt-d2003-main-reclassify-startwithakm",
            mean  ="m2-mixt-d2003-main-reclassify-startwithmean",
            va    ="m2-mixt-d2003-main-reclassify-startwithva",
            prank = "m2-mixt-d2003-main-reclassify-startwithprank",
            mratio= "m2-mixt-d2003-main-reclassify-startwithratiomover",
            resid = "m2-mixt-d2003-main-reclassify-residuals",
            dynanmic = "m4-mixt-d2003-dirty-reclassification")
  for (nn in names(ll)) {
    if (!res.exists(ll[[nn]])) next;
    res = res.load(ll[[nn]])
    rr = ldply(res,function(v) data.frame(v$vdec$stats))
    rr$start= nn
    rr_all = rbind(rr_all,rr)
  }
  res.save("m2-mixt-d2003-main-reclassify-allresults",rr_all)
}

server.static.mixt.estimate.fit.bs <- function() {

  load(sprintf("%s/data-tmp/data-static.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]
  sim = list(jdata=jdata,sdata=sdata)
  rm(sdata,jdata)

  # get clsuters
  grps  = res.load("m2-mixt-y2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  # get main estimate
  res_main = res.load("m2-mixt-y2003-main-fixb")

  # limit to firms with 5 workers
  fids = sim$sdata[,.N,f1][N>=5,f1]
  sim$jdata = sim$jdata[f1 %in% fids][f2 %in% fids]
  sim$sdata = sim$sdata[f1 %in% fids]

  # impute data with BLM
  jdata.sim = m2.mixt.impute.movers(sim$jdata,res_main$model)[,list(y1=y1_imp,j1,j2,y2=y2_imp,f1,f2)]
  sdata.sim = m2.mixt.impute.stayers(sim$sdata,res_main$model)
  sdata.sim[,y1:=y1_imp][,y2:=y2_imp]
  sim2 = list(sdata=sdata.sim,jdata=jdata.sim)

  m2.bs(sim2)
  m2.bs(sim)

}

# Werun the whole BLM procedure but only using residuals
server.static.mixt.estimate.robustess.residuals <- function() {

  sim = server.static.data()
  grps  = res.load("m2-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  sdata = sim$sdata
  jdata = sim$jdata

  # create some variables
  # get industry for all firms
  firm_info = unique(sdata[,list(f1,ind1)])
  setkey(firm_info,f1)
  setkey(sdata,f2)
  sdata[,ind2:= firm_info[sdata,ind1]]
  sdata[,educ_f := factor(educ)]
  jdata[,educ_f := factor(educ)]
  sdata[,age:=2002 - birthyear]
  jdata[,age:=2002 - birthyear]

  # compute a wage regression or wages
  adata = rbind(sdata[,list(wid,age,ind=ind1,educ_f,t=1,y=y1)],sdata[,list(wid,age,ind=ind2,educ_f,t=2,y=y2)])

  fit = lm(y ~ I(age)*educ_f + I(age^2)*educ_f + I(age^3)*educ_f + ind*educ_f*I(t==1) ,adata,na.action = na.exclude)

  # predict wages for sdata, jdata
  sdata[, y1t := y1]
  sdata[, y1  := y1t - predict(fit,newdata=sdata[,list(age,ind=ind1,educ_f,t=1)])]
  sdata[, y2t := y2]
  sdata[, y2  := y2t - predict(fit,newdata=sdata[,list(age,ind=ind2,educ_f,t=2)])]
  jdata[, y1t := y1]
  jdata[, y1  := y1t - predict(fit,newdata=jdata[,list(age,ind=ind1,educ_f,t=1)])]
  jdata[, y2t := y2]
  jdata[, y2  := y2t - predict(fit,newdata=jdata[,list(age,ind=ind2,educ_f,t=2)])]


  sim = list(jdata=jdata,sdata=sdata[move==FALSE])
  rm(sdata,jdata)

  # we start with clustering
  set.seed(65542134)
  ms    = grouping.getMeasures(sim,"ecdf",Nw=20,y_var = "y1")
  grps  = grouping.classify.once(ms,k = 10,nstart = 1000,iter.max = 200,step=100)
  sim   = grouping.append(sim,grps$best_cluster)

  ctrl      = em.control(nplot=1000,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=local_opts$estimation.mixture$est_rep,
                         est_nbest=local_opts$estimation.mixture$est_nbest,
                         sdata_subsample=0.1,
                         maxiter = local_opts$estimation.mixture$maxiter)

  cl = makeCluster(local_opts$cpu_count)
  clusterEvalQ(cl,require(blmrep))

  res_mixt = m2.mixt.estimate.all(sim,nk=6,ctrl,cl)
  stopCluster(cl)

  res.save("m2-mixt-d2003-main-residuals",res_mixt)
}

#' This is for the second part of the main results for static
server.static.analysis.meaneffects <- function() {
  res_main     = res.load("m2-mixt-d2003-main-fixb") # we do it across bootstraps
  res_bs       = res.load("m2-mixt-d2003-bootstrap")

  NNs          = res_main$model$NNs
  NNs          = round(NNs*1e6/sum(NNs))

  rr = data.frame()
  for (nn in names(res_bs$mixt_all)) {
    res = res_bs$mixt_all[[nn]]
    rt = m2.mixt.meaneffect(res$model)
    rr = rbind(rr,rt)
    flog.info("done with %s",nn)
  }

  rr  = data.table(rr)
  rrm = melt(rr,id.vars = "type")
  rrs = rrm[, list(m1=mean(value), q0 = quantile(value,0.025), q1 = quantile(value,0.975),sd = sd(value)) , list(variable,type)]

  #  --------- attach the main result --------- #
  # finally attach the main results
  rt = m2.mixt.meaneffect(res_main$model)
  rrs[type=="pk0",main := as.numeric(rt[type=="pk0",1:7])]
  rrs[type=="pku",main := as.numeric(rt[type=="pku",1:7])]

  #rt = m2.mixt.meaneffect(res_main_hl$model)
  rrs = rrs[,list(variable, type, main, m1,q0,q1,sd)]

  # ---- do in differences ----- #
  rrd  = rrm[,value:= value[type=="pku"] - value[type=="pk0"] , list(variable)]
  rrds = rrd[, list(m1=mean(value), q0 = quantile(value,0.025), q1 = quantile(value,0.975),sd = sd(value)) , list(variable)]
  rrds[,main := as.numeric(rt[type=="pku",1:7]) - as.numeric(rt[type=="pk0",1:7]) ]
  rrds = rrds[,list(variable, main, m1,q0,q1,sd)]

  res.save("m2-mixt-d2003-meffects",list(level=rrs,diff=rrds))
}

# estimates the mixture of mixture model
server.static.mixture.mixtofmixt <- function() {

  sim   = server.static.data()
  grps  = res.load("m2-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  # get main estimate
  res_main = res.load("m2-mixt-d2003-main-fixb")
  model_mixt = res_main$model
  model_mixt$S1 = 0.9*model_mixt$S1 + 0.1*mean(model_mixt$S1)
  model_mixt$S2 = 0.9*model_mixt$S2 + 0.1*mean(model_mixt$S2)
  model_np   = m2.mixt.np.new.from.ns(model_mixt,3)

  ctrl   = em.control(tol=1e-6,fixb=F,dprior=1.001,ncat=1,posterior_reg=1e-7,sd_floor=1e-6,
                      maxiter = local_opts$estimation.mixture$maxiter)

  # estimate movers
  res_np = m2.mixt.np.movers.estimate(sim$jdata,model_np,ctrl)
  res_np$model$pk0 = res_np$model$pk0[1,,] # not using covariates.

  # estimate stayers
  sim$sdata[,sample := rank(runif(.N))/.N<=ctrl$sdata_subsample,j1]
  res_np = m2.mixt.np.stayers.estimate(sim$sdata[sample==1],res_np$model,ctrl)
  # rr = m2.mixt.np.movers.residuals(res_np$model,TRUE,100)

  stayer_share = sim$sdata[,.N]/(sim$sdata[,.N]+sim$jdata[,.N])
  vdec = m2.mixt.np.vdec(res_np$model,nsim = 1e6,stayer_share = stayer_share)
  res_np$vdec = vdec
  res.save("m2-mixt-d2003-main-mixtofmixt",res_np)
}


# ========== PROBABILISTIC APPROACH =====

server.static.proba.results <- function() {

  sim   = server.static.data()
  grps  = res.load("m2-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  # ----------- Prepare starting classification --------
  # ru AKM on jdata to extract firm FE
  firm.fe = m2.fe.firms(sim$jdata)

  # classify according to AKM
  require(Ckmeans.1d.dp)
  clusters   = Ckmeans.1d.dp(firm.fe$fe$psi, 10)
  grp_akm        = clusters$cluster
  names(grp_akm) = firm.fe$fe$f1

  sim   = server.static.data()
  sim   = grouping.append(sim,grp_akm,drop=T)

  # Setting up the options for BLM
  ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=local_opts$estimation.mixture$est_rep,
                         est_nbest=local_opts$estimation.mixture$est_nbest,
                         sdata_subsample=0.1,
                         maxiter = local_opts$estimation.mixture$maxiter)

  cl = makeCluster(local_opts$cpu_count)
  clusterEvalQ(cl,require(blmrep))
  res_mixt_akm = m2.mixt.estimate.all(sim,nk=6,ctrl,cl)
  res.save("m2-akmgrp-d2003-fixb",res_mixt_akm)

  # using the same sample, cluster and estimate
  ms           = grouping.getMeasures(sim,"ecdf",Nw=20,y_var = "y1")
  grps_kmean   = grouping.classify.once(ms,k = 10,nstart = 1000,iter.max = 200,step=100)
  sim          = grouping.append(sim,grps_kmean$best_cluster)
  res_mixt_blm = m2.mixt.estimate.all(sim,nk=6,ctrl,cl)
  stopCluster(cl)

  res.save("m2-akmgrp-d2003-fixb",list(res_mixt_akm=res_mixt_akm,res_mixt_blm=res_mixt_blm,grps_kmean=grps_kmean$best_cluster,grp_akm=grp_akm))

  # ----------- Gibbs iteration statrting with AKM --------
  res_all = res.load("m2-akmgrp-d2003-fixb")
  sim      = server.static.data()
  iter_res = res_all$res_mixt_akm
  sim      = grouping.append(sim,res_all$grp_akm,drop=T)

  iter_pl  = m2.proba.getlcasspr(sim,10)$class_pr
  iter_grp = res_all$grp_akm
  # loop re-calssifcation/estimation
  rr.all  = data.frame()
  rr.all2 = data.frame()
  for (i in 1:local_opts$estimation.probabilistic$maxiter) {

    # run gibbs 200 firms
    iter_gibbs = m2.proba.gibbs(sim,iter_res$model,iter_pl,nfirms_to_update = 1,grp_start = iter_grp,maxiter = local_opts$estimation.probabilistic$gibbs_nfirm)
    iter_grp   = iter_gibbs$last_grp

    # apply latest classification
    sim        = grouping.append(sim,iter_gibbs$last_grp)
    iter_pl    = m2.proba.getlcasspr(sim,10)$class_pr

    # estimate BLM using nodes
    cl = makeCluster(local_opts$cpu_count)
    clusterEvalQ(cl,require(blmrep))
    iter_res = m2.mixt.estimate.all(sim,nk=6,ctrl,cl)
    stopCluster(cl)

    rr.all = rbind(rr.all,data.frame(iter_res$vdec$stats))
    iter_gibbs$rr$outloop=i
    rr.all2 = rbind(rr.all2,iter_gibbs$rr)
  }

  gibbs.all = list()
  gibbs.all$akm$vdec = rr.all
  gibbs.all$akm$liks = rr.all2

  # ----------- Gibbs iteration statrting with BLM --------
  res_all = res.load("m2-akmgrp-d2003-fixb")
  sim      = server.static.data()
  iter_res = res_all$res_mixt_blm
  sim      = grouping.append(sim,res_all$grps_kmean,drop=T)

  iter_pl  = m2.proba.getlcasspr(sim,10)$class_pr
  iter_grp = res_all$grps_kmean
  # loop re-calssifcation/estimation
  rr.all  = data.frame()
  rr.all2 = data.frame()
  for (i in 1:local_opts$estimation.probabilistic$maxiter) {

    # run gibbs 200 firms
    iter_gibbs = m2.proba.gibbs(sim,iter_res$model,iter_pl,nfirms_to_update = 1,grp_start = iter_grp,maxiter = local_opts$estimation.probabilistic$gibbs_nfirm)
    iter_grp   = iter_gibbs$last_grp

    # apply latest classification
    sim        = grouping.append(sim,iter_gibbs$last_grp)
    iter_pl    = m2.proba.getlcasspr(sim,10)$class_pr

    # estimate BLM using nodes
    cl = makeCluster(local_opts$cpu_count)
    clusterEvalQ(cl,require(blmrep))
    iter_res = m2.mixt.estimate.all(sim,nk=6,ctrl,cl)
    stopCluster(cl)

    rr.all = rbind(rr.all,data.frame(iter_res$vdec$stats))
    iter_gibbs$rr$outloop=i
    rr.all2 = rbind(rr.all2,iter_gibbs$rr)
  }

  gibbs.all$blm$vdec = rr.all
  gibbs.all$blm$liks = rr.all2

  res.save("m2-proba-gibbs-d2003-res",gibbs.all)
}


# ========== FIXED EFFECT MODEL ===========
server.fe.trace <- function() {
  sim = server.static.data()
  res = m2.trace.blmhyb(sim,use_con = TRUE,nm_list = local_opts$trace$nm_list)
  res.save("m2-blmhybrid",res)
}


# ========== SHIMER-SMITH with OTJ  ===========
server.shimersmith.results <- function() {
  p <- blmrep:::initp( b=0.3, c=0, sz=1, nx=6, ny = 10)
  p$pf = function(x,y,z=0,p)  ( 0.5*x^p$rho +  0.5*y^p$rho )^(1/p$rho) + p$ay # haggerdorn law manovski
  p$ay  = 0.7 #0.5
  p$rho = -3 # 2.5 #-1.5 # 2.5 # -1.5

  r = list()
  p$rho = -3 # 2.5 #-1.5 # 2.5 # -1.5
  r$pam_6x10 = shimersmith.simulateAndEstimate(p,est_rep = local_opts$estimation.mixture$est_rep,
                                               est_nbest = local_opts$estimation.mixture$est_nbest,
                                               maxiter = local_opts$estimation.mixture$maxiter)
  p$rho = 3  # 2.5 #-1.5 # 2.5 # -1.5
  r$nam_6x10 = shimersmith.simulateAndEstimate(p,est_rep = local_opts$estimation.mixture$est_rep,
                                               est_nbest = local_opts$estimation.mixture$est_nbest,
                                               maxiter = local_opts$estimation.mixture$maxiter)
  res.save("m2-shimersmith-mc",r)
}


