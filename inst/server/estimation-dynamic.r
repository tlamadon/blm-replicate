
server.dynamic.data <- function(remove_movers_from_sdata=T) {
  load(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))
  sdata <- sdata[,x:=1][,y1_bu:=y1]
  if (remove_movers_from_sdata) {
    sdata=sdata[move==FALSE]
  }
  sim = list(jdata=jdata,sdata=sdata)
  rm(sdata,jdata)
  return(sim)
}

# ==== DYNAMIC DATA =====

server.dynamic.d2003.computeclusters <- function() {

  sim = server.dynamic.data()

  set.seed(954938257)
  ms    = grouping.getMeasures(sim,"ecdf",Nw=20,y_var = "y2")
  grps  = grouping.classify.once(ms,k = 10,nstart = 10000,iter.max = 200,step=100)
  sim   = grouping.append(sim,grps$best_cluster)
  res.save("m4-mixt-d2003-groups",grps)
}

server.dynamic.d2003.clustering.stats <- function() {

  sim   = server.dynamic.data()
  grps  = res.load("m4-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  mstats = sim$jdata[,list(m1=mean(y1),sd1=sd(y1),
                           m2=mean(y2),sd2=sd(y2),
                           m1=mean(y3),sd1=sd(y3),
                           m1=mean(y4),sd1=sd(y4),.N),list(j1,j2)]
  cstats = sim$sdata[,list(m1=mean(y1),sd1=sd(y1),
                           m2=mean(y2),sd2=sd(y2),
                           m1=mean(y3),sd1=sd(y3),
                           m1=mean(y4),sd1=sd(y4),.N),list(j1)]

  sdata = sim$sdata
  setkey(sdata,j1)
  sdata[,age:=2002 - birthyear]
  sdata[,asize := .N,f1]
  sdata[,fmw := mean(y2),f1]

  gstats = rBind(
    sdata[, data.frame(get.stats.clusters(.SD,ydep="y2")), j1],
    data.frame(list(get.stats.clusters(sdata,ydep="y2"),j1=0)))

  res.save("m4-stats-d2003",list(mstats=mstats,cstats=cstats,gstats = gstats))
}

# === MINI MODEL 2003 ====

# used in try/catch
err_fun <- function(e) {catf("error in boot strap rep %i!\n",i);print(e)}



server.static.rho.analysis <- function(){

  load(sprintf("%s/data-tmp/tmp-2003-static",local_opts$wdir))
  clus = unique(sdata[,list(fid=f1,clus=j1)])
  load(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))
  sdata = cluster.append.data(sdata,clus)
  jdata = cluster.append.data(jdata,clus)

  res_rhos = model.mini4.getvar.stayers.unc2.opt(sdata)

  # get slices for each rho
  ff <- function(rhos) model.mini4.getvar.stayers.unc2(sdata,rhos[1],rhos[2],rhos[3],as.numeric(weights))$fit1

  rr = data.frame()
  rsup = seq(0.01,0.99,l=100)
  vals = sapply(rsup,function(x) ff(c(x,res_rhos$r4,res_rhos$rt)))
  rr = rbind(rr,data.frame(x=rsup,y=vals,name="rho1",best=res_rhos$r1))
  vals = sapply(rsup,function(x) ff(c(res_rhos$r1,x,res_rhos$rt)))
  rr = rbind(rr,data.frame(x=rsup,y=vals,name="rho4",best=res_rhos$r4))
  vals = sapply(rsup,function(x) ff(c(res_rhos$r1,res_rhos$r4,x)))
  rr = rbind(rr,data.frame(x=rsup,y=vals,name="rho32",best=res_rhos$rt))

  #ggplot(rr,aes(x=x,y=y)) + geom_line() + facet_grid(~name) + theme_bw() + geom_vline(aes(xintercept=best))
  res.save("rho-analysis",rr)
}

#' estimate the mini-model once, then resimulate from it, then re-cluster
#' and re-estimate. repeat theis many times
server.dynamic.mini.estimate <-function() {

  sim   = server.dynamic.data()
  grps  = res.load("m4-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  set.seed(12345)

  # estimate
  res_rhos    = m4.mini.getvar.stayers.unc2.opt(sim$sdata)
  model_mini  = m4.mini.estimate(sim$jdata,sim$sdata,res_rhos$r1,res_rhos$r4,method="prof")
  model_lin   = m4.mini.estimate(sim$jdata,sim$sdata,res_rhos$r1,res_rhos$r4,method="linear")

  res.save("m4-mini-d2003")
}

#' estimate the mini-model with different cluster sizes
server.dynamic.mini.estimate.ksize <- function() {

  load(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))

  rr_models = list()
  rr_vardecs =list()
  for (k in 3:20) {
    jdata =  cluster.append.data(jdata,res[[paste(k)]])
    sdata =  cluster.append.data(sdata,res[[paste(k)]])
    res_rhos_bis        = model.mini4.getvar.stayers.unc2.opt(sdata)
    model_mini_bis      = model.mini4.estimate(jdata,sdata,res_rhos_bis$r1,res_rhos_bis$r4,fixb=T)
    sdata.sim.sim       = model.mini4.impute.stayers(model_mini_bis,sdata)
    proj_unc_mini_bis1  = lin.proja(sdata.sim.sim,"y2_imp","k_imp","j1");
    proj_unc_mini_bis1$between_var_explained = sdata[, list(mean(y2),.N),j1][,wtd.var(V1,N)]/sdata[, list(mean(y2),.N),f1][,wtd.var(V1,N)]*100
    rr_models[[paste(k)]] = model_mini_bis
    rr_vardecs[[paste(k)]] = proj_unc_mini_bis1
  }

  archive.put(mini_ksize = rr_models ,m="mini model for different k sizes",file = arch_dyn)
  archive.put(mini_ksize_vardec=rr_vardecs,m="all variance decomposition for each cluster size",file = arch_dyn)
}

server.dynamic.cf <- function(){

  load(sprintf("%s/data-tmp/tmp-2003-static",local_opts$wdir))
  clus = unique(sdata[,list(fid=f1,clus=j1)])
  load(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)
  clus2 = clus$clus
  names(clus2)=clus$fid
  sim = grouping.append(sim,clus2)

  mstats = sim$jdata[,list(m1=mean(y1),sd1=sd(y1),
                       m2=mean(y2),sd2=sd(y2),
                       m1=mean(y3),sd1=sd(y3),
                       m1=mean(y4),sd1=sd(y4),.N),list(j1,j2)]
  cstats = sim$sdata[,list(m1=mean(y1),sd1=sd(y1),
                       m2=mean(y2),sd2=sd(y2),
                       m1=mean(y3),sd1=sd(y3),
                       m1=mean(y4),sd1=sd(y4),.N),list(j1)]

  res_rhos = m4.mini.getvar.stayers.unc2.opt(sim$sdata)
  model_mini      = m4.mini.estimate(sim$jdata,sim$sdata,res_rhos$r1,res_rhos$r4,method="prof")
  m4.mini.plotw(model_mini)

  sdata.sim       = m4.mini.impute.stayers(model_mini,sim$sdata)
  vdec_minimodel  = lin.proja(sdata.sim,"y2_imp","k_imp","j1");

  sdata.sim[,move:=FALSE]
  jdata.sim       = model.mini4.impute.movers(model_mini,jdata)
  data.sim=rbind(sdata.sim,jdata.sim)
  reg=list()
  reg[["ana_y2_move"]]      = coef((data.sim[,lm(y2_imp~factor(j1)+move+k_imp)]))
  reg[["ana_y2_move_i"]]    = coef((data.sim[,lm(y2_imp~factor(j1)*k_imp+move)]))
  reg[["ana_y2_k1k2"]]      = coef((jdata.sim[,lm(y2_imp~factor(j1)+factor(j2)+k_imp)]))
  reg[["ana_y3_k1k2"]]      = coef((jdata.sim[,lm(y3_imp~factor(j1)+factor(j2)+k_imp)]))
  reg[["ana_y3_k1k2y2"]]    = coef((jdata.sim[,lm(y3_imp~factor(j1)+factor(j2)+k_imp+y2_imp)]))
  reg[["ana_y3_k1k2y2_i"]]  = coef((jdata.sim[,lm(y3_imp~factor(j1)+factor(j2)*k_imp+y2_imp)]))

}

#' estimate the mini-model once, then resimulate from it, then re-cluster
#' and re-estimate. repeat theis many times
server.dynamic.mini.estimate.bootstrap <-function() {

  load(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)

  set.seed(12345)

  ms    = grouping.getMeasures(sim,"ecdf",Nw=20,y_var = "y2")
  grps  = grouping.classify(ms,ksupp = 10,stop = TRUE,nstart = 1000,verbose =10,iter.max = 200)
  sim   = grouping.append(sim,grps$best_cluster)

  # estimate
  res_rhos = m4.mini.getvar.stayers.unc2.opt(sim$sdata)
  model0   = m4.mini.estimate(sim$jdata,sim$sdata,res_rhos$r1,res_rhos$r4,method="prof")
  # m4.mini.plotw(model0)
  sdata.sim   = m4.mini.impute.stayers(model0,sim$sdata)
  vdec2 = lin.proja(sdata.sim,"y2_imp","k_imp","j1")

  dclus0 = sim$sdata[,list(j1=j1[1]),f1]
  clus0  = dclus0$j1
  names(clus0) = dclus0$f1

  rrr=data.frame()
  rr_mix = list()
  rr_mini = list()
  for (i in 1:100) {
    sim   = grouping.append(sim,clus0,sort = F)
    sim$sdata = m4.mini.impute.stayers(model0,sim$sdata)
    sim$jdata = m4.mini.impute.movers(model0,sim$jdata)
    sim$sdata[,c('y1','y2','y3','y4'):=list(y1_imp,y2_imp,y3_imp,y4_imp)]
    sim$jdata[,c('y1','y2','y3','y4'):=list(y1_imp,y2_imp,y3_imp,y4_imp)]

    ms    = grouping.getMeasures(sim,"ecdf",Nw=10,y_var = "y2")
    grps  = grouping.classify(ms,ksupp = 10,stop = TRUE,nstart = 200,verbose =10,iter.max = 200)
    sim   = grouping.append(sim,grps$best_cluster)

    res_rhos         = m4.mini.getvar.stayers.unc2.opt(sim$sdata)
    mini_model_bis   = m4.mini.estimate(sim$jdata,sim$sdata,res_rhos$r1,res_rhos$r4,method="prof")

    sdata.sim   = m4.mini.impute.stayers(mini_model_bis,sim$sdata[sample.int(.N,1e5)])
    vdec2 = lin.proja(sdata.sim,"y2_imp","k_imp","j1")
    rr_mini[[paste(i)]] = mini_model_bis
    rr2 = data.frame(vdec2$stats)
    rr2$model="mini"
    rr2$i= i

    catf("don with rep %i\n",i)
    rrr = rbind(rrr,rr2)
  }

  gg = data.table(ldply(rr_mini,function(x) x$B1/x$B1[2]))
  mg = melt(gg,id.vars = ".id")[,list(m=mean(value),sd=sd(value)),variable]

  # need to boostrap they actual figure :-/
  gg = data.table(ldply(rr_mini,function(x) m4.mini.plotw(x,getvals=T)))
  gg = gg[,list(value=mean(value),ql=quantile(value,0.05),qh=quantile(value,0.95)),list(l,k,variable)]
  ggplot(gg,aes(x=l,y=value,color=factor(k))) + geom_line() +
    geom_point() + geom_errorbar(aes(ymin=ql,ymax=qh),width=0.2)+
    theme_bw() + facet_wrap(~variable,nrow = 2,scales="free") + coord_cartesian(ylim=c(9.3,11.5))

  # compute the covariance
  gg = data.table(ldply(rr_mini,function(x) cov.wt(cbind(x$B1,x$Em),x$Ns)$cov[2,1]))
  gg[,list(m=mean(V1),q0=quantile(V1,0.025),q1=quantile(V1,0.975))]
  cov.wt(cbind(model0$B1,model0$Em),model0$Ns)$cov[2,1]

  I = 2:10
  gg = data.table(ldply(rr_mini,function(x) cov.wt(cbind(x$B1[I],x$Em[I]),x$Ns[I])$cov[2,1]))
  gg[,list(m=mean(V1),q0=quantile(V1,0.025),q1=quantile(V1,0.975))]
  cov.wt(cbind(model0$B1[I],model0$Em[I]),model0$Ns[I])$cov[2,1]

}

#' Explores the sensitivity to rho
server.dynamic.mini.estimate.explore_rho <- function(){
  load("L:\\Tibo\\qtrdata\\tmp-2003-static.dat")
  clus = unique(sdata[,list(fid=f1,clus=j1)])
  load(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))
  sdata = cluster.append.data(sdata,clus)
  jdata = cluster.append.data(jdata,clus)

  res_rhos             = model.mini4.getvar.stayers.unc2.opt(sdata)
  mini_model_true      = model.mini4.estimate(jdata,sdata,res_rhos$r1,res_rhos$r4,fixb=T)

  # generate values at which we are going to evaluate the model
  pw = 4
  rr = data.frame(r1 = seq( (0.3*res_rhos$r1)^pw,0.95^pw,l=20)^(1/pw),r4 = seq( (0.3*res_rhos$r4)^pw,0.95^pw,l=20)^(1/pw))
  rrr = data.frame()
  for (i in 1:nrow(rr)) {
    model_mini      = model.mini4.estimate(jdata,sdata,rr$r1[i],rr$r4[i],fixb=T)
    sdata.sim       = model.mini4.impute.stayers(model_mini,sdata)
    vdec_minimodel  = lin.proja(sdata.sim,"y2_imp","k_imp","j1");
    rtmp = data.frame(vdec_minimodel$stats)
    rtmp$r1=rr$r1[i]
    rtmp$r4=rr$r4[i]
    rtmp$neg_b = sum(model_mini$B1<0)
    rrr = rbind(rrr,rtmp)
  }
}

# ======= MIXTURE MODEL 2003 =============

server.dynamic.mixture.d2003.estimate <- function() {
  sim = server.dynamic.data()
  grps = res.load("m4-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster,drop = T)

  set.seed(54140598)
  cl = makeCluster(local_opts$cpu_count)
  clusterEvalQ(cl,require(blmrep))

  ctrl      = em.control(est_rho=c(FALSE,TRUE,FALSE,FALSE,TRUE,FALSE),
                         nplot=10000,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=local_opts$estimation.mixture$est_rep,
                         est_nbest=local_opts$estimation.mixture$est_nbest,
                         sdata_subsample=0.1,
                         maxiter = local_opts$estimation.mixture$maxiter)
  res_mixt  = m4.mixt.estimate.all(sim,nk=6,ctrl,cl)
  stopCluster(cl)
  res.save("m4-mixt-d2003-main",res_mixt)
}

server.dynamic.mixture.d2003.estimate.model_iteration <- function() {
  sim = server.dynamic.data()
  grps = res.load("m4-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster,drop = T)

  ctrl      = em.control(nplot=2000,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=local_opts$estimation.mixture$est_rep,
                         est_nbest=local_opts$estimation.mixture$est_nbest,
                         sdata_subsample=0.1,
                         maxiter = local_opts$estimation.mixture$maxiter)
  cl = makeCluster(local_opts$cpu_count)
  clusterEvalQ(cl,require(blmrep))
  res   = m4.mixt.estimate.reclassify(sim,maxiter=10,ctrl=ctrl,cl=cl)
  stopCluster(cl)

  res.save("m4-mixt-d2003-reclassify",res)
}

server.dynamic.mixture.d2003.fit <- function() {
  load(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))
  sdata[,y1_bu := y1]
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)

  # get the groups
  grps = res.load("m4-mixt-2003data-groups")
  sim   = grouping.append(sim,grps$best_cluster)
  sim$sdata[,list(y1,y2,y3,y4)]

  # get main estimate
  res_main = res.load("m4-mixt-d2003-main")

  # impute data on movers
  jdata.sim = m4.mixt.impute.movers(sim$jdata,res_main$model)

  vmin = quantile(jdata.sim$y2,0.01)
  vmax = quantile(jdata.sim$y2,0.99)

  jdata.sim[, j1b := ceiling(j1/2)]
  jdata.sim[, j2b := ceiling(j2/2)]

  dd_m = jdata.sim[, {
    d1 = density(y2,bw="sj",from=vmin,to=vmax,n=30)
    d2 = density(y2_imp,bw="sj",from=vmin,to=vmax,n=30)
    data.frame(x=d1$x,y1=d1$y,y2=d2$y)
  }, list(j1b=sprintf("%i-%i",2*j1b-1,2*j1b),j2b=sprintf("%i-%i",2*j2b-1,2*j2b)) ]

  # growth
  vmin = quantile(jdata.sim[,y3-y2],0.01)
  vmax = quantile(jdata.sim[,y3-y2],0.99)

  dd_m_g = jdata.sim[, {
    d1 = density(y3-y2,bw="sj",from=vmin,to=vmax,n=30)
    d2 = density(y3_imp-y2_imp,bw="sj",from=vmin,to=vmax,n=30)
    data.frame(x=d1$x,y1=d1$y,y2=d2$y)
  }, list(j1b=sprintf("%i-%i",2*j1b-1,2*j1b),j2b=sprintf("%i-%i",2*j2b-1,2*j2b)) ]

  rr_m = m4.movers.checkfit(jdata.sim,res_main$model$B12,res_main$model$B43)

  sdata.sim = m4.mixt.impute.stayers(sim$sdata,res_main$model)
  rr_s      = m4.stayers.checkfit(sdata.sim,res_main$model$B12,res_main$model$B43)

  vmin = quantile(sdata.sim$y1,0.001)
  vmax = quantile(sdata.sim$y1,0.999)

  dd_s = sdata.sim[, {
    d1 = density(y2,from=vmin,to=vmax,n=30)
    d2 = density(y2_imp,from=vmin,to=vmax,n=30)
    data.frame(x=d1$x,y1=d1$y,y2=d2$y)
  }, list(j1) ]

  # saving the distribubtion of growth
  vmin = quantile(sdata.sim[,y2-y1],0.001)
  vmax = quantile(sdata.sim[,y2-y1],0.999)

  dd_s_g = sdata.sim[, {
    d1 = density(y2-y1,from=vmin,to=vmax,n=30)
    d2 = density(y2_imp-y1_imp,from=vmin,to=vmax,n=30)
    data.frame(x=d1$x,y1=d1$y,y2=d2$y)
  }, list(j1) ]

  rs    = rkiv0.start("m4-mixt-d2003-main-fit",info="fit of the model")
  res.save(rs,list(dd_s=dd_s,dd_m=dd_m,rr_s=rr_s,rr_m=rr_m,dd_s_g=dd_s_g,dd_m_g=dd_m_g))
}

#' this estimates the model with different values of K to check
#' robustness. We recluster every time
server.dynamic.mixture.estimate.robust.nf <- function() {

  sim = server.dynamic.data()
  grps = res.load("m4-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster,drop = T)

  set.seed(87954352)
  ctrl      = em.control(nplot=2000,tol=1e-7,dprior=1.001,fixb=TRUE,
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
      ms    = grouping.getMeasures(sim,"ecdf",Nw=20,y_var = "y2")
      grps  = grouping.classify.once(ms,k = nf_size,nstart = 1000,iter.max = 200)
      sim   = grouping.append(sim,grps$best_cluster)
      res_mixt = m4.mixt.estimate.all(sim,nk=6,ctrl,cl=cl)
      rr_mixt[[paste("nf",nf_size,sep="-")]] = res_mixt
    })
  }
  stopCluster(cl)

  res.save("m4-mixt-d2003-change-nf",rr_mixt)
}


#' this estimates the model with different values of K to check
#' robustness. We recluster every time
server.dynamic.mixture.estimate.robust.nk <- function() {
  sim = server.dynamic.data()
  grps = res.load("m4-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster,drop = T)

  set.seed(87954352)
  ctrl      = em.control(nplot=2000,tol=1e-7,dprior=1.001,fixb=TRUE,
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
      res_mixt = m4.mixt.estimate.all(sim,nk=nk_size,ctrl=ctrl,cl=cl)
      res_mixt$second_stage_reps_all=NA
      rr_mixt[[paste("nk",nk_size,sep="-")]] = res_mixt
      save(rr_mixt,file="tmp-m4-chg-nk.dat")
    })}
  stopCluster(cl)

  res.save("m4-mixt-d2003-change-nf",rr_mixt)
}

server.dynamic.mixture.estimate.robust.different_rho <- function() {
  sim = server.dynamic.data()
  grps = res.load("m4-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster,drop = T)

  set.seed(87954352)
  ctrl      = em.control(nplot=2000,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=local_opts$estimation.mixture$est_rep,
                         est_nbest=local_opts$estimation.mixture$est_nbest,
                         sdata_subsample=0.1,
                         maxiter = local_opts$estimation.mixture$maxiter)

  cl = makeCluster(local_opts$cpu_count)
  clusterEvalQ(cl,require(blmrep))

  res_rhos   = m4.mini.getvar.stayers.unc2.opt(sim$sdata,diff=TRUE)

  res_mixt = m4.mixt.estimate.all(sim,nk=6,ctrl=ctrl,cl=cl)
  stopCluster(cl)

  res.save("m4-mixt-d2003-rho-check",res_mixt)
}

server.dynamic.mixture.d2003.boostrap <- function(){

  sim   = server.dynamic.data()
  grps  = res.load("m4-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster,drop = T)

  set.seed(879543542)
  ctrl      = em.control(est_rho=c(FALSE,TRUE,FALSE,FALSE,TRUE,FALSE),
                         nplot=1000,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-5,posterior_reg=1e-8,
                         est_rep=local_opts$estimation.mixture$est_rep,
                         est_nbest=local_opts$estimation.mixture$est_nbest,
                         sdata_subsample=0.1,
                         maxiter = local_opts$estimation.mixture$maxiter)

  # get main estimate
  res_main = res.load("m4-mixt-d2003-main")

  cl = makeCluster(local_opts$cpu_count)
  clusterEvalQ(cl,require(blmrep))

  nfirms=0
  rr = data.frame()
  rr_mixt_all = list()
  for (rep in 1:local_opts$bootstrap_nreps) {

    sdata.sim = m4.mixt.impute.stayers(sim$sdata,res_main$model)
    jdata.sim = m4.mixt.impute.movers(sim$jdata,res_main$model)
    sdata.sim = sdata.sim[,y1:=y1_imp][,y2:=y2_imp][,y3:=y3_imp][,y4:=y4_imp][,jt:=j1][,y1_bu:=y1]
    jdata.sim = jdata.sim[,y1:=y1_imp][,y2:=y2_imp][,y3:=y3_imp][,y4:=y4_imp]
    sim.sp = list(jdata=jdata.sim,sdata=sdata.sim)
    rm(jdata.sim,sdata.sim)

    # re-cluster
    ms              = grouping.getMeasures(sim.sp,"ecdf",Nw=20,y_var = "y2")
    grps            = grouping.classify(ms,ksupp=10, stop=FALSE, nstart= 500, verbose=10, iter.max=200)
    sim.sp          = grouping.append(sim.sp,grps$best_cluster,drop = T)

    # get the rhos
    sim.sp$sdata[,y1:=y1_bu]
    res_rhos        = m4.mini.getvar.stayers.unc2.opt(sim.sp$sdata)

    # run mixture
    res_mixt = m4.mixt.estimate.all(sim.sp,nk=6,ctrl,cl=cl)
    rt = data.frame(res_mixt$vdec$stats)
    rt$name="mixt";rt$rep=rep;rt$nfirms=nfirms
    rt$r1=res_mixt$model$B12
    rt$r4=res_mixt$model$B43
    rr = rbind(rr,rt)

    res_mixt$clusi_cor = sim.sp$sdata[,cor(jt,j1)]
    res_mixt$second_stage_reps_all = NULL
    rr_mixt_all[[paste(rep)]] = res_mixt
    rm(sim.sp)

    print(rr)
  }

  stopCluster(cl)
  res.save("m4-mixt-d2003-bootstrap",list(vdecs=rr,mixt_all=rr_mixt_all))
}

# ============== computing counter factuals ===============

#' This is for the second part of the main results for dynamic
server.dynamic.analysis.meaneffects <- function() {
  res_main    = res.load("m4-mixt-d2003-main")
  res_bs      = res.load("m4-mixt-d2003-bootstrap")

  rr = data.frame()
  for (nn in names(res_bs$mixt_all)) {
    res = res_bs$mixt_all[[nn]]
    rt = m4.mixt.meaneffect(res$model)
    rr = rbind(rr,rt)
    flog.info("done with %s",nn)
  }

  rr  = data.table(rr)
  rrm = melt(rr,id.vars = "type")
  rrs = rrm[, list(m1=mean(value), q0 = quantile(value,0.025), q1 = quantile(value,0.975),sd = sd(value)) , list(variable,type)]

  # -----------  finally attach the main results --------- #
  rt = m4.mixt.meaneffect(res_main$model)
  #rt = m4.mixt.meaneffect(res_main_hl$model)
  rrs[type=="pk0",main := as.numeric(rt[type=="pk0",1:7])]
  rrs[type=="pku",main := as.numeric(rt[type=="pku",1:7])]
  rrs = rrs[,list(variable, type, main, m1,q0,q1)]

  # ---- do in differences ----- #
  rrd  = rrm[,value:= value[type=="pku"] - value[type=="pk0"] , list(variable)]
  rrds = rrd[, list(m1=mean(value), q0 = quantile(value,0.025), q1 = quantile(value,0.975),sd = sd(value)) , list(variable)]
  rrds[,main := as.numeric(rt[type=="pku",1:7]) - as.numeric(rt[type=="pk0",1:7]) ]
  rrds = rrds[,list(variable, main, m1,q0,q1,sd)]

  res.save("m4-mixt-d2003-meffects",list(level=rrs,diff=rrds))
}

#' This is the endogenous mobility results
server.dynamic.analysis.endogenousmobility <- function() {

  # get the main result
  res_main      = res.load("m4-mixt-d2003-main")
  model         = res_main$model
  dd = compute.mob.matrix(res_main$model)

  # get boostrap replications
  res_bs      = res.load("m4-mixt-d2003-bootstrap")
  rr = ldply(res_bs$mixt_all, function(x) {
    compute.mob.matrix(x$model)},.progress = "text")
  rr = data.table(rr)

  drr = rr[.id>100,list(m=mean(value),q0=quantile(value,0.025),q1=quantile(value,0.975),sd=sd(value)),list(j1,j2,cond)]
  res.save("m4-cf-mobility",list(stats=drr,all=rr,main=dd))
}

server.dynamic.analysis.stateDependence <- function() {

  res_main    = res.load("m4-mixt-d2003-main")
  res_bs      = res.load("m4-mixt-d2003-bootstrap")
  model = res_main$model

  # combine everything
  rr = analysis.dynamic.dec.bis(model,1e6)

  # run the bootstrap
  rr.bs = ldply(res_bs$mixt_all,function(m) analysis.dynamic.dec.bis(m$model,1e6),.progress = "text")
  rrm = data.table(melt(rr.bs,id=".id"))
  rrs = rrm[, list(m1=mean(value,na.rm=T), q0 = quantile(value,0.025,na.rm=T), q1 = quantile(value,0.975,na.rm=T),sd = sd(value,na.rm=T)) , list(variable)]

  rrs[,main:=as.numeric(rr)]
  rrs = rrs[,list(variable,main,m1,sd,q0,q1)]

  res.save("m4-cf-state-dependence",rrs)
}







