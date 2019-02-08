
server.dynamic.data <- function(remove_movers_from_sdata=T) {
  load(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]
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




server.dynamic.mixture.d2003.estimate.old <- function() {

  load(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))
  sdata[,y1_bu := y1]
  sdata = sdata[move==FALSE]
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)

  # get the groups
  grps = res.load("m4-mixt-2003data-groups")
  sim   = grouping.append(sim,grps$best_cluster)
  sim$sdata[,list(y1,y2,y3,y4)]

  set.seed(54140598)
  rs      = rkiv0.start("m4-mixt-d2003-main")
  rs$info = "dynamic mixture on 2003 data"

  ctrl      = em.control(est_rho=c(FALSE,TRUE,FALSE,FALSE,TRUE,FALSE),
                       nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                       sd_floor=1e-7,posterior_reg=1e-8,
                       est_rep=50,est_nbest=10,sdata_subsample=0.1)
  res_mixt  = m4.mixt.estimate.all(sim,nk=6,ctrl,cl)
  res.save(rs,res_mixt)

  # ------ decomposition on best likelihood ----- #
  # get main estimate
  res_main  = res.load("m4-mixt-d2003-main")
  rrm.order = res_main$second_stage_reps
  name = rrm.order[order(-lik_mixt)][1,i]
  res  = res_main$second_stage_reps_all[[name]]
  sim$sdata[,S1 := rank(runif(.N))/.N<=ctrl$sdata_subsample,j1]
  sim$sdata[,y1:=y1_bu]
  sim$sdata[,S1 := rank(runif(.N))/.N<=ctrl$sdata_subsample,j1]
  res = m4.mixt.rhoext.stayers(sim$jdata, sim$sdata[S1==1], res$model, em.control(ctrl,cstr_type="akm",cstr_val=0.05,fixb=FALSE))
  sim$sdata[,y1:=y1_bu]
  res = m4.mixt.rhoext.stayers(sim$jdata, sim$sdata[S1==1], res$model, em.control(ctrl,cstr_type="none"))

  proj = m4.mixt.vdec(res$model,1e6,sim$sdata[,.N]/(sim$sdata[,.N]+sim$jdata[,.N]),"y2")
  res$vdec = proj

  rs    = rkiv0.start("m4-mixt-d2003-main-highestlik",
                      info="dynamic 2003 main estimate, vdec on highest likelihood of movers")
  res.save(rs,res)



  # ---- all best 10 liks ------ #
  res_mixt = res.load("m4-mixt-d2003-main")
  rrm.order = res_mixt$second_stage_reps
  rrs = list()
  for (ii in 1:10) {
    name = rrm.order[order(-lik_mixt)][ii,i]
    res  = res_mixt$second_stage_reps_all[[name]]
    sim$sdata[,y1:=y1_bu]
    sim$sdata[,S1 := rank(runif(.N))/.N<=ctrl$sdata_subsample,j1]
    res = m4.mixt.rhoext.stayers(sim$jdata, sim$sdata[move==FALSE][S1==1], res$model, em.control(ctrl,cstr_type="akm",cstr_val=0.05,fixb=FALSE))
    sim$sdata[,y1:=y1_bu]
    res = m4.mixt.rhoext.stayers(sim$jdata, sim$sdata[move==FALSE][S1==1], res$model, em.control(ctrl,cstr_type="none"))

    # simulate movers/stayers, and combine
    NNm = res$model$NNm
    NNs = res$model$NNs/ctrl$sdata_subsample
    NNm[!is.finite(NNm)]=0
    NNs[!is.finite(NNs)]=0
    share_s  = sum(NNs)/(sum(NNm) + sum(NNs))
    share_m  = sum(NNm)/(sum(NNm) + sum(NNs))

    NNs = round(NNs*ctrl$vdec_sim_size*share_s/sum(NNs))
    NNm = round(NNm*ctrl$vdec_sim_size*share_m/sum(NNm))

    # we simulate from the model both movers and stayers
    sdata.sim = m4.mixt.simulate.stayers(res$model,NNs)
    jdata.sim = m4.mixt.simulate.movers(res$model,NNm)
    sdata.sim = rbind(sdata.sim[,list(j1,k,y2)],jdata.sim[,list(j1,k,y2)])
    proj_unc  = lin.proj(sdata.sim,"y2","k","j1");

    res$proj = proj_unc
    res$tau=NULL
    rrs[[name]]=res
  }

  rs      = rkiv0.start("m4-mixt-d2003-main-stayer-reps")
  rs$info = "variance decompositions for 10 best starting values"
  res.save(rs,rrs)

  # ---- MINI MODEL ---- #
  # running mini model unconstrained ad with B=1 on same data
  set.seed(690923542)
  rs         = rkiv0.start("m4-mini-d2003-profandlinear")
  rs$info    = "dynamic mini model on pooled data, using profiling and with B=1"
  sim$sdata[,y1:=y1_bu]
  res_rhos   = m4.mini.getvar.stayers.unc2.opt(sim$sdata)

  sim$sdata[,y1:=y1_bu]
  model_mini = m4.mini.estimate(sim$jdata,sim$sdata,res_rhos$r1,res_rhos$r4,method="prof")
  sim$sdata[,y1:=y1_bu]
  model_minilin      = m4.mini.estimate(sim$jdata,sim$sdata,res_rhos$r1,res_rhos$r4,method="linear")
  res.save(rs,list(mini=model_mini,mini_lin=model_minilin))
}

server.dynamic.mixture.d2003.estimate <- function() {
  sim = server.dynamic.data()
  grps = res.load("m4-mixt-2003data-groups")
  sim   = grouping.append(sim,grps$best_cluster,drop = T)

  set.seed(54140598)
  cl = makeCluster(15)
  clusterEvalQ(cl,require(blmrep))

  ctrl      = em.control(est_rho=c(FALSE,TRUE,FALSE,FALSE,TRUE,FALSE),
                         nplot=10000,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=50,est_nbest=10,sdata_subsample=0.1)
  res_mixt  = m4.mixt.estimate.all(sim,nk=6,ctrl,cl)
  stopCluster(cl)
  res.save("m4-mixt-d2003-main",res_mixt)

}

server.dynamic.mixture.d2003.estimate.reclassify <- function() {

  load(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))
  sdata[,y1_bu := y1]
  sdata = sdata[move==FALSE]
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)

  # get the groups
  grps = res.load("m4-mixt-2003data-groups")
  sim   = grouping.append(sim,grps$best_cluster)
  sim$sdata[,list(y1,y2,y3,y4)]

  res = m4.mixt.estimate.reclassify(sim,maxiter=10)
  rs      = rkiv0.start("m4-mixt-d2003-dirty-reclassification")
  rs$info = "running reclassificaiton on dynamic model"
  res.save(rs,res)

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


# here we retain 10% of the stayers in classification
server.static.mixture.d2003.estimate.holdout <- function() {

  load(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]
  sdata = sdata[move==FALSE]
  sim = list(jdata=jdata,sdata=sdata)
  rm(sdata,jdata)

  rs    = rkiv0.start("m4-mixt-d2003-main-holdout",
                      info="dynamic 2003 main estimate, holding stayers in clustering")

  # create holdout
  sim$sdata[,HO := rank(runif(.N))/.N <= 0.1 ,f1]

  ms    = grouping.getMeasures(list(sdata=sim$sdata[HO==FALSE],jdata=sim$jdata),"ecdf",Nw=20,y_var = "y2")
  grps  = grouping.classify.once(ms,k = 10,nstart = 1000,iter.max = 200,step=100)
  sim   = grouping.append(sim,grps$best_cluster)

  ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=34,est_nbest=10,sdata_subsample=1,sdata_subredraw=FALSE)
  cl = makeCluster(17)
  clusterEvalQ(cl,require(blm))

  sim$sdata[,sample:=HO]
  res_mixt = m4.mixt.estimate.all(sim,nk=6,ctrl,cl)
  stopCluster(cl)

  res.save(rs,res_mixt)

  # ---- split stayers ----- #
  # try to include half the stayers in clustering
  load(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]

  # split movers in 2
  jdata[,splitdata := rank(runif(.N)) - 0.5 +0.1*(runif(1)-0.5) <= .N/2,f1]
  mwids = jdata[splitdata==TRUE,unique(wid)]
  sdata = sdata[!(wid %in% mwids)]
  jdata = jdata[wid %in% mwids]
  sim = list(jdata=jdata,sdata=sdata)
  rm(sdata,jdata)

  rs    = rkiv0.start("m4-mixt-d2003-main-holdout2",
                      info="dynamic 2003 main estimate, using 50% of movers in clustering")

  ms    = grouping.getMeasures(sim,"ecdf",Nw=20,y_var = "y2")
  grps  = grouping.classify.once(ms,k = 10,nstart = 1000,iter.max = 200,step=100)
  sim   = grouping.append(sim,grps$best_cluster)

  ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=50,est_nbest=10,sdata_subsample=0.1)
  cl = makeCluster(17)
  clusterEvalQ(cl,require(blm))
  res_mixt = m4.mixt.estimate.all(list(sdata=sim$sdata[move==FALSE],jdata=sim$jdata),nk=6,ctrl,cl)
  stopCluster(cl)

  res.save(rs,res_mixt)





}

server.dynamic.mixture.estimate.2003data.split <- function() {

  load(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))
  sdata[,y1_bu := y1]
  sdata = sdata[move==FALSE] # do not include movers in sdata
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)

  # get the groups
  grps = res.load("m4-mixt-2003data-groups")
  sim   = grouping.append(sim,grps$best_cluster)
  sim$sdata[,list(y1,y2,y3,y4)]

  ctrl      = em.control(est_rho=c(FALSE,TRUE,FALSE,FALSE,TRUE,FALSE),
                         nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=50,est_nbest=10,sdata_subsample=0.1)

  no_cores = 17
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl,require(blm))

  # --- SPLITTING SAMPLEs ----- #
  set.seed(13928498)
  rs    = rkiv0.start("m4-mixt-d2003-splits",dependson="m4-mixt-d2003-groups",
                      info="dynamic 2003 main estimate")
  sim$sdata[,splitdata := rank(runif(.N)) - 0.5 +0.1*(runif(1)-0.5) <= .N/2,f1] # splitting stayers
  sim$jdata[,splitdata := rank(runif(.N)) - 0.5 +0.1*(runif(1)-0.5) <= .N/2,f1] # splitting movers

  sim.sp          = list(sdata= sim$sdata[splitdata==TRUE],jdata = sim$jdata[splitdata==TRUE])
  ms              = grouping.getMeasures(sim.sp,"ecdf",Nw=20,y_var = "y2")
  grps            = grouping.classify.once(ms,k10, nstart= 500, verbose=10, iter.max=200,step=50)
  sim.sp          = grouping.append(sim.sp,grps$best_cluster,drop = T)
  res_mixt_split1 = m4.mixt.estimate.all(sim.sp,nk=6,ctrl,cl=cl)
  model_mini1     = m4.mini.estimate(sim.sp$jdata,sim.sp$sdata,res_mixt_split1$model$B12,res_mixt_split1$model$B43,method="prof")
  model_minilin1  = m4.mini.estimate(sim.sp$jdata,sim.sp$sdata,res_mixt_split1$model$B12,res_mixt_split1$model$B43,method="linear")

  sim.sp          = list(sdata= sim$sdata[splitdata==FALSE],jdata = sim$jdata[splitdata==FALSE])
  ms              = grouping.getMeasures(sim.sp,"ecdf",Nw=20,y_var = "y2")
  grps            = grouping.classify.once(ms,k10, nstart= 500, verbose=10, iter.max=200,step=50)
  sim.sp          = grouping.append(sim.sp,grps$best_cluster,drop=T)
  res_mixt_split2 = res_mixt = m4.mixt.estimate.all(sim.sp,nk=6,ctrl,cl=cl)
  model_mini2     = m4.mini.estimate(sim.sp$jdata,sim.sp$sdata,res_mixt_split2$model$B12,res_mixt_split2$model$B43,method="prof")
  model_minilin2  = m4.mini.estimate(sim.sp$jdata,sim.sp$sdata,res_mixt_split2$model$B12,res_mixt_split2$model$B43,method="linear")

  stopCluster(cl)

  res.save(rs,list(
    split1=list(mixt=res_mixt_split1,mini=model_mini1,minilin=model_minilin1),
    split2=list(mixt=res_mixt_split2,mini=model_mini2,minilin=model_minilin2)))

}


#' this estimates the model with different values of K to check
#' robustness. We recluster every time
server.dynamic.mixture.estimate.robust.nf <- function() {

  load(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))
  sdata[,y1_bu := y1]
  sdata = sdata[move==FALSE]
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)

  set.seed(87954352)
  rs    = rkiv0.start("m4-mixt-y2003-changeK",
                      info="dynamic 2003 estimation with different number of clusters")
  ctrl      = em.control(nplot=2000,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=45,est_nbest=10,sdata_subsample=0.1)

  cl = makeCluster(15)
  clusterEvalQ(cl,require(blm))

  rr_mixt = list()
  for (nf_size in c(3:15,17,20)) {
    tryCatch({
      ms    = grouping.getMeasures(sim,"ecdf",Nw=20,y_var = "y2")
      grps  = grouping.classify.once(ms,k = nf_size,nstart = 1000,iter.max = 200)
      sim   = grouping.append(sim,grps$best_cluster)
      res_mixt = m4.mixt.estimate.all(sim,nk=6,ctrl,cl=cl)
      rr_mixt[[paste("nf",nf_size,sep="-")]] = res_mixt
    })}
  stopCluster(cl)

  res.save(rs,rr_mixt)
}


#' this estimates the model with different values of K to check
#' robustness. We recluster every time
server.dynamic.mixture.estimate.robust.nk <- function() {

  load(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))
  sdata[,y1_bu := y1]
  sdata = sdata[move==FALSE]
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)

  # get the groups
  grps = res.load("m4-mixt-2003data-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  set.seed(87954352)
  rs    = rkiv0.start("m4-mixt-d2003-change-nk",
                      info="dynamic 2003 estimation with different number of worker types")
  ctrl      = em.control(nplot=2000,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=45,est_nbest=10,sdata_subsample=0.1)

  cl = makeCluster(15)
  clusterEvalQ(cl,require(blm))

  rr_mixt = list()
  rr_mixt = list()
  for (nk_size in c(3:5,7:9)) {
    tryCatch({
      res_mixt = m4.mixt.estimate.all(sim,nk=nk_size,ctrl=ctrl,cl=cl)
      res_mixt$second_stage_reps_all=NA
      rr_mixt[[paste("nk",nk_size,sep="-")]] = res_mixt
      save(rr_mixt,file="tmp-m4-chg-nk.dat")
    })}
  stopClsuter(cl)
  res.save(rs,rr_mixt)

}

server.dynamic.mixture.estimate.robust.different.rho <- function() {

  load(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))
  sdata[,y1_bu := y1]
  sdata = sdata[move==FALSE]
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)

  # get the groups
  grps = res.load("m4-mixt-2003data-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  ctrl      = em.control(nplot=2000,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=45,est_nbest=10,sdata_subsample=0.1,rho_in_diff=TRUE)

  cl = makeCluster(15)
  clusterEvalQ(cl,require(blm))
  res_mixt = m4.mixt.estimate.all(sim,nk=6,ctrl=ctrl,cl=cl)
  stopCluster(cl)

  rs    = rkiv0.start("m4-mixt-d2003-rho-check",
                      info="estimating rho in differences")
  res.save(rs,res_mixt_rho)

  res_rhos   = m4.mini.getvar.stayers.unc2.opt(sim$sdata,diff=TRUE)
  sim$sdata[,y1:=y1_bu]
  model_mini = m4.mini.estimate(sim$jdata,sim$sdata,res_rhos$r1,res_rhos$r4,method="prof")

  ddres = data.table(data=res_rhos$dep,model=res_rhos$fitted,w = res_rhos$weights)
  ggplot(ddres,aes(x=data,y=model,size=w)) +
    geom_point() + geom_abline(linetype=2) + theme_bw()

  # check the fit of the covariance matrix
  res_mixt_rho = res.load("m4-mixt-d2003-rho-check")

  sdata2 = m4.mixt.impute.stayers(sim$sdata,res_mixt_rho$model)
  sdata2[,var(cbind(y1,y2,y3,y4))]
  sdata2[,var(cbind(y1_imp,y2_imp,y3_imp,y4_imp))]
  sdata2[,var(cbind(y2-y1,y3-y2,y4-y3))]
  sdata2[,var(cbind(y2_imp-y1_imp,y3_imp-y2_imp,y4_imp-y3_imp))]
}

server.dynamic.mixture.d2003.boostrap <- function(){

  load(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))
  sdata[,y1_bu := y1]
  sdata=sdata[move==FALSE]
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)

  # get the groups
  grps = res.load("m4-mixt-2003data-groups")
  sim   = grouping.append(sim,grps$best_cluster)
  sim$sdata[,list(y1,y2,y3,y4)]

  # get main estimate
  res_main = res.load("m4-mixt-d2003-main")

  ctrl      = em.control(est_rho=c(FALSE,TRUE,FALSE,FALSE,TRUE,FALSE),
                         nplot=1000,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-5,posterior_reg=1e-8,
                         est_rep=45,est_nbest=10,sdata_subsample=0.1)

  no_cores = 15
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl,require(blm))

  # --- SPLITTING SAMPLEs ----- #
  nfirms=0

  rr = data.frame()
  rr_mixt_all = list()

  load("tmpbootstrap-dyn-bis-lap3.dat")

  reps = setdiff(1:100,as.integer(names(rr_mixt_all)))
  for (rep in reps) {
    sdata.sim = m4.mixt.impute.stayers(sim$sdata,res_main$model)
    jdata.sim = m4.mixt.impute.movers(sim$jdata,res_main$model)
    sdata.sim = sdata.sim[,y1:=y1_imp][,y2:=y2_imp][,y3:=y3_imp][,y4:=y4_imp][,jt:=j1][,y1_bu:=y1]
    jdata.sim = jdata.sim[,y1:=y1_imp][,y2:=y2_imp][,y3:=y3_imp][,y4:=y4_imp]
    sim.sp = list(jdata=jdata.sim,sdata=sdata.sim)
    rm(jdata.sim,sdata.sim)

    sim.sp$sdata[,y1:=y1_bu]
    res_rhos        = m4.mini.getvar.stayers.unc2.opt(sim.sp$sdata)
    model_mini1     = m4.mini.estimate(sim.sp$jdata,sim.sp$sdata,res_rhos$r1,res_rhos$r4,method="prof")
    model_minilin1  = m4.mini.estimate(sim.sp$jdata,sim.sp$sdata,res_rhos$r1,res_rhos$r4,method="linear")

    rt = data.frame(model_mini1$vdec$stats)
    rt$name="mini-truecluster";rt$rep=rep;rt$nfirms=nfirms
    rt$r1=res_rhos$r1
    rt$r4=res_rhos$r4
    rr = rbind(rr,rt)
    rt = data.frame(model_minilin1$vdec$stats)
    rt$name="linear-truecluster";rt$rep=rep;rt$nfirms=nfirms
    rt$r1=res_rhos$r1
    rt$r4=res_rhos$r4
    rr = rbind(rr,rt)

    # re-cluster
    ms              = grouping.getMeasures(sim.sp,"ecdf",Nw=20,y_var = "y2")
    grps            = grouping.classify(ms,ksupp=10, stop=FALSE, nstart= 500, verbose=10, iter.max=200)
    sim.sp          = grouping.append(sim.sp,grps$best_cluster,drop = T)

    sim.sp$sdata[,y1:=y1_bu]
    res_rhos        = m4.mini.getvar.stayers.unc2.opt(sim.sp$sdata)
    model_mini1     = m4.mini.estimate(sim.sp$jdata,sim.sp$sdata,res_rhos$r1,res_rhos$r4,method="prof")
    model_minilin1  = m4.mini.estimate(sim.sp$jdata,sim.sp$sdata,res_rhos$r1,res_rhos$r4,method="linear")

    rt = data.frame(model_mini1$vdec$stats)
    rt$name="mini";rt$rep=rep;rt$nfirms=nfirms
    rt$r1=res_rhos$r1
    rt$r4=res_rhos$r4
    rr = rbind(rr,rt)

    rt = data.frame(model_minilin1$vdec$stats)
    rt$name="linear";rt$rep=rep;rt$nfirms=nfirms
    rt$r1=res_rhos$r1
    rt$r4=res_rhos$r4
    rr = rbind(rr,rt)

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

    save(rr,rr_mixt_all,file = "tmpbootstrap-dyn-bis-lap3.dat")

    print(rr)
  }

  stopCluster(cl)

  load("tmpbootstrap-dyn-bis-lap3.dat")
  rs      = rkiv0.start("m4-mixt-d2003-bootstrap-leg2")
  rs$info = "parametric bootstrap of dynamic mixture on 2003 data, reps 101 to 200"
  res.save(rs,list(vdecs=rr,mixt_all=rr_mixt_all))


}

server.dynamic.mixture.d2003.boostrap.resample <- function(){

  load(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))
  sdata[,y1_bu := y1]
  sdata=sdata[move==FALSE]
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)

  # get the groups
  grps = res.load("m4-mixt-2003data-groups")
  sim   = grouping.append(sim,grps$best_cluster)
  sim$sdata[,list(y1,y2,y3,y4)]

  # get main estimate
  res_main = res.load("m4-mixt-d2003-main")

  ctrl      = em.control(est_rho=c(FALSE,TRUE,FALSE,FALSE,TRUE,FALSE),
                         nplot=1000,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-5,posterior_reg=1e-8,
                         est_rep=45,est_nbest=10,sdata_subsample=0.1)

  no_cores = 15
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl,require(blm))

  # --- SPLITTING SAMPLEs ----- #
  nfirms=0

  rr = data.frame()
  rr_mixt_all = list()
  rep_start=1

  #load("tmpbootstrap-dyn-bis.dat")
  #rep_start = max(rr$rep)+1
  Ntot = sim$jdata[,.N] + sim$sdata[,.N]

  for (rep in rep_start:100) {
    # we resample
    wids = sample(c(sim$jdata$wid,sim$sdata$wid),Ntot,replace=T)

    setkey(sim$sdata,wid)
    I = (wids %in% c(sim$sdata$wid))
    sdata.sim = sim$sdata[wids[I]]
    setkey(sim$jdata,wid)
    I = (wids %in% c(sim$jdata$wid))
    jdata.sim = sim$jdata[wids[I]]

    sdata.sim = m4.mixt.impute.stayers(sdata.sim,res_main$model)
    jdata.sim = m4.mixt.impute.movers(jdata.sim,res_main$model)
    sdata.sim = sdata.sim[,y1:=y1_imp][,y2:=y2_imp][,y3:=y3_imp][,y4:=y4_imp][,jt:=j1][,y1_bu:=y1]
    jdata.sim = jdata.sim[,y1:=y1_imp][,y2:=y2_imp][,y3:=y3_imp][,y4:=y4_imp]
    sim.sp = list(jdata=jdata.sim,sdata=sdata.sim)
    rm(jdata.sim,sdata.sim)

    sim.sp$sdata[,y1:=y1_bu]
    res_rhos        = m4.mini.getvar.stayers.unc2.opt(sim.sp$sdata)
    model_mini1     = m4.mini.estimate(sim.sp$jdata,sim.sp$sdata,res_rhos$r1,res_rhos$r4,method="prof")
    model_minilin1  = m4.mini.estimate(sim.sp$jdata,sim.sp$sdata,res_rhos$r1,res_rhos$r4,method="linear")

    rt = data.frame(model_mini1$vdec$stats)
    rt$name="mini-truecluster";rt$rep=rep;rt$nfirms=nfirms
    rt$r1=res_rhos$r1
    rt$r4=res_rhos$r4
    rr = rbind(rr,rt)
    rt = data.frame(model_minilin1$vdec$stats)
    rt$name="linear-truecluster";rt$rep=rep;rt$nfirms=nfirms
    rt$r1=res_rhos$r1
    rt$r4=res_rhos$r4
    rr = rbind(rr,rt)

    # re-cluster
    ms              = grouping.getMeasures(sim.sp,"ecdf",Nw=20,y_var = "y2")
    grps            = grouping.classify.once(ms,k=10, nstart= 1000, verbose=10, iter.max=200)
    sim.sp          = grouping.append(sim.sp,grps$best_cluster,drop = T)

    sim.sp$sdata[,y1:=y1_bu]
    res_rhos        = m4.mini.getvar.stayers.unc2.opt(sim.sp$sdata)
    model_mini1     = m4.mini.estimate(sim.sp$jdata,sim.sp$sdata,res_rhos$r1,res_rhos$r4,method="prof")
    model_minilin1  = m4.mini.estimate(sim.sp$jdata,sim.sp$sdata,res_rhos$r1,res_rhos$r4,method="linear")

    rt = data.frame(model_mini1$vdec$stats)
    rt$name="mini";rt$rep=rep;rt$nfirms=nfirms
    rt$r1=res_rhos$r1
    rt$r4=res_rhos$r4
    rr = rbind(rr,rt)

    rt = data.frame(model_minilin1$vdec$stats)
    rt$name="linear";rt$rep=rep;rt$nfirms=nfirms
    rt$r1=res_rhos$r1
    rt$r4=res_rhos$r4
    rr = rbind(rr,rt)

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

    save(rr,rr_mixt_all,file = "tmpbootstrap-dyn-resample.dat")

    print(rr)
  }

  stopCluster(cl)

  load("tmpbootstrap-dyn-resample.dat")
  rs      = rkiv0.start("m4-mixt-d2003-bootstrap-resample")
  rs$info = "parametric bootstrap of dynamic mixture on 2003 data"
  res.save(rs,list(vdecs=rr,mixt_all=rr_mixt_all))


}


server.dynamic.d2003.changeK <- function() {

  load(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))
  sdata[,y1_bu := y1]
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)

  rs    = rkiv0.start("m4-mixt-d2003-changeK",
                      info="dynamic 2003 estimation with different number of clusters")
  ctrl      = em.control(est_rho=c(FALSE,TRUE,FALSE,FALSE,TRUE,FALSE),
                         nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=45,est_nbest=10,sdata_subsample=0.1)

  no_cores = 15
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl,require(blm))

  rr_mixt = list()
  for (nf_size in c(6,8,12,15)) {
      ms    = grouping.getMeasures(sim,"ecdf",Nw=20,y_var = "y2")
      grps  = grouping.classify(ms,ksupp = nf_size,stop = TRUE,nstart = 500,verbose =10,iter.max = 200)
      sim   = grouping.append(sim,grps$best_cluster)
      res_mixt       = m4.mixt.estimate.all(sim,nk=6,ctrl,cl=cl)
      model_mini     = m4.mini.estimate(sim$jdata,sim$sdata,res_rhos$r1,res_rhos$r4,method="prof")
      model_minilin  = m4.mini.estimate(sim$jdata,sim$sdata,res_rhos$r1,res_rhos$r4,method="linear")
      rr_mixt[[paste("nf",nf_size,sep="-")]] = list(mixt = res_mixt, mini=model_mini,mini_lin=model_mini_lin)
  }
  res.save(rs,rr_mixt)
  stopCluster(cl)
}




server.dynamic.mixture.estimate.2003data.localmax <- function() {

  load(sprintf("%s/data-tmp/tmp-2003-dynamic.dat",local_opts$wdir))
  sdata[,y1_bu := y1]
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)

  # get the groups
  grps = res.load("m4-mixt-2003data-groups")
  sim   = grouping.append(sim,grps$best_cluster)
  sim$sdata[,list(y1,y2,y3,y4)]

  # get best estimates
  res_mixt = res.load("m4-mixt-2003data-main-para-bis")

  # try data perturbation
  ctrl      = em.control(est_rho=c(FALSE,TRUE,FALSE,FALSE,TRUE,FALSE),
                         nplot=500,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,ncat=1,
                         est_rep=50,est_nbest=10,sdata_subsample=0.1)
  ctrl$pert_type="data"
  ctrl$pert_beta=10
  res2 = m4.mixt.rhoext.movers.pert(sim$jdata,sim$sdata,res_mixt$model,ctrl = ctrl)

  ctrl$pert_type="sa"
  ctrl$pert_beta=0.01
  res2 = m4.mixt.rhoext.movers.pert(sim$jdata,sim$sdata,res_mixt$model,ctrl = ctrl)
  for (i in 1:20) {
    ctrl$pert_beta=ctrl$pert_beta^0.7
    res2 = m4.mixt.rhoext.movers.pert(sim$jdata,sim$sdata,res2$model,ctrl = ctrl)
  }


  # try to strart from para with error to see if we get higher
  res = res.load("m4-mixt-2003data-main-para")
  ctrl      = em.control(est_rho=c(FALSE,TRUE,FALSE,FALSE,TRUE,FALSE),
                         nplot=500,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,ncat=1,
                         est_rep=50,est_nbest=10,sdata_subsample=0.1)
  res2 = m4.mixt.rhoext.movers(sim$jdata,sim$sdata,res$model,ctrl = ctrl)


  # start from the high variance result, amd re-estimate
  res = res.load("m4-mixt-2003data-main")
  ctrl      = em.control(est_rho=c(FALSE,TRUE,FALSE,FALSE,TRUE,FALSE),
                         nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,ncat=1,
                         est_rep=50,est_nbest=10,sdata_subsample=0.1)
  res2 = m4.mixt.rhoext.movers(sim$jdata,sim$sdata,res$model,ctrl = ctrl)

  ctrl$pert_type="sa"
  ctrl$pert_beta=0.4
  res2 = m4.mixt.rhoext.movers.pert(sim$jdata[sample.int(.N,ceiling(.N/2))],sim$sdata,res$model,ctrl = ctrl)
  for (i in 1:20) {
    ctrl$pert_beta=ctrl$pert_beta^0.5
    flog.info("beta=%f",ctrl$pert_beta)
    res2 = m4.mixt.rhoext.movers.pert(sim$jdata,sim$sdata,res2$model,ctrl = ctrl)
  }
  res2 = m4.mixt.rhoext.movers(sim$jdata,sim$sdata,res2$model,ctrl = ctrl)

  # compute decomposition
  sim$sdata[,S1 := rank(runif(.N))/.N<=ctrl$sdata_subsample,j1]
  sim$sdata[,y1:=y1_bu]
  res2 = m4.mixt.rhoext.stayers(sim$jdata, sim$sdata[move==FALSE][S1==1], res2$model, em.control(ctrl,cstr_type="akm",cstr_val=0.05))
  sim$sdata[,y1:=y1_bu]
  res2 = m4.mixt.rhoext.stayers(sim$jdata, sim$sdata[move==FALSE][S1==1], res2$model, em.control(ctrl,cstr_type="none"))
  sim$sdata[,y1:=y1_bu]
  sdatae.sim = m4.mixt.impute.stayers( sim$sdata[move==FALSE][S1==1],res$model)
  proj_unc   = blm:::lin.proj(sdatae.sim[k_imp!=5],"y2_imp","k_imp","j1");

  model.connectiveness3(res2$model,T)

}


server.dynamic.mini.explore <- function() {

  fmeas <- function(dt) {
      res = lm(y2~y1,dt);
      v = c(mean(dt$y2),sd(dt$y2),coef(res)[1],coef(res)[2]);
      v[is.na(v)]=0
      data.frame(m=c("m1","m2","m3","m4"),value=v)
  }

  fmeas <- function(dt) {
    res = lm(y2~y1,dt);
    v = c(mean(dt$y2),sd(dt$y2),cov(dt$y1,dt$y2));
    v[is.na(v)]=0
    data.frame(m=c("m1","m2","m3"),value=v)
  }

  fmeas <- function(dt) {
    v = c(mean(dt$y1 -res_rhos$r1 * dt$y2 ),sd(dt$y1 -res_rhos$r1 * dt$y2 ));
    v[is.na(v)]=0
    data.frame(m=c("m1","m2"),value=v)
  }

  meas2 = cluster.getMeasures(sdata,jdata,measure="user_dyn",user_fn=fmeas,normalize = TRUE)

  source("../../aclus/aclus/R/firm-clustering.R")

  grp2 = firm.group(meas2,ksupp=10)
  jdata[,j1:= as.integer(grp2$best_cluster[[f1]]),f1]
  jdata[,j2:= as.integer(grp2$best_cluster[[f2]]),f2]
  sdata[,j1:=NULL]
  sdata[,j1:= as.integer(grp2$best_cluster[[f1]]),f1]

  acast(jdata[,.N,list(j1,j2)],j1~j2)

  res_rhos = model.mini4.getvar.stayers.unc2.opt(sdata)
  model_mini      = model.mini4.estimate(jdata,sdata,res_rhos$r1,res_rhos$r4,fixb=T)
  sdata.sim       = model.mini4.impute.stayers(model_mini,sdata)
  proj_unc_mini   = lin.proja(sdata.sim,"y2_imp","k_imp","j1");

  rbind(model_mini$B2,model_mini$A1)

  # check the fit
  sdata[,as.list(coef(lm(y2~y1))),j1]
  sdata.sim[,as.list(coef(lm(y2~y1))),j1]


}







