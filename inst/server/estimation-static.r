server.static.data <- function(remove_movers_from_sdata=T) {
  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
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
  grps  = rkiv0.load("m2-mixt-y2003-groups")
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
}

#' estimate the mini-model once, then resimulate from it, then re-cluster
#' and re-estimate. repeat theis many times
server.static.mini.estimate.bootstrap <-function() {

  sim = server.static.data()
  grps  = res.load("m2-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  # get the main mini-model estimate
  model0 = res.load("m2-mini-prof")
  sdata.sim   = m2.mini.impute.stayers(model0,sim$sdata)
  vdec2 = lin.proja(sdata.sim,"y1_imp","k_imp","j1")

  # save original clusters
  dclus0 = unique(sim$sdata[,list(f1,j1)])
  clus0  = dclus0$j1
  names(clus0) = dclus0$f1

  # set the seed
  set.seed(12345) # archive should save that too

  rrr=data.frame()
  rr_mix = list()
  rr_mini = list()
  for (i in 1:100) {
    sim   = grouping.append(sim,clus0,sort = F)
    sim$sdata = m2.mini.impute.stayers(model0,sim$sdata)
    sim$jdata = m2.mini.impute.movers(model0,sim$jdata)
    sim$sdata[,c('y1','y2'):=list(y1_imp,y2_imp)]
    sim$jdata[,c('y1','y2'):=list(y1_imp,y2_imp)]

    ms    = grouping.getMeasures(sim,"ecdf",Nw=10,y_var = "y1")
    grps  = grouping.classify(ms,ksupp = 10,stop = TRUE,nstart = 300,verbose =10,iter.max = 200)
    sim   = grouping.append(sim,grps$best_cluster)

    model1   = m2.mini.estimate(sim$jdata,sim$sdata,method="prof",norm=4)
    # m2.mini.plotw(model1)

    sdata.sim   = m2.mini.impute.stayers(model1,sim$sdata)
    vdec2 = lin.proja(sdata.sim,"y1_imp","k_imp","j1")
    rr_mini[[paste(i)]] = model1
    rr2 = data.frame(vdec2$stats)
    rr2$model="mini"
    rr2$i= i

    catf("don with rep %i\n",i)
    rrr = rbind(rrr,rr2)
  }

  gg = data.table(ldply(rr_mini,function(x) x$B1/x$B1[2]))
  mg = melt(gg,id.vars = ".id")[,list(m=mean(value),sd=sd(value)),variable]

  # need to boostrap they actual figure :-/
  gg = data.table(ldply(rr_mini,function(x) m2.mini.plotw(x,getvals=T)))
  gg = gg[,list(value=mean(value),ql=quantile(value,0.025),qh=quantile(value,0.975)),list(l,k,variable)]
  ggplot(gg,aes(x=l,y=value,color=factor(k))) + geom_line() +
    geom_point() + geom_errorbar(aes(ymin=ql,ymax=qh),width=0.2)+
    theme_bw() + facet_wrap(~variable,nrow = 2,scales="free")

  # compute the covariance
  gg = data.table(ldply(rr_mini,function(x) cov.wt(cbind(x$B1,x$Em),x$Ns)$cov[2,1]))
  gg[,list(m=mean(V1),q0=quantile(V1,0.025),q1=quantile(V1,0.975))]
  cov.wt(cbind(model0$B1,model0$Em),model0$Ns)$cov[2,1]

  I = 2:10
  gg = data.table(ldply(rr_mini,function(x) cov.wt(cbind(x$B1[I],x$Em[I]),x$Ns[I])$cov[2,1]))
  gg[,list(m=mean(V1),q0=quantile(V1,0.025),q1=quantile(V1,0.975))]

  archive.put(mini_model_bs=rr_mini,m="main mini boostrapped, 2003 data",file = arch_static)
  archive.put(mini_model_bs_vdec=rrr,m="variance decomposition for main mini boostrapped, 2003",file = arch_static)

  # compute decompositions
  data.table(melt(rrr,id.vars = c("model","i")))[,list(m=mean(value),q0=quantile(value,0.025),q1=quantile(value,0.975)),variable]


}

server.static.mini.estimate.mainx <- function(){
  arch_static = "res-2003-static.dat"

  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,ageg:= (age<=30) + 2*((age >= 31)&(age<=50)) + 3*(age>50)]
  sdata[,x := educ + 3*(ageg-1)]
  sdata[, x := as.integer(x)]

  cstats = sdata[,list(m1=mean(y1),sd1=sd(y1),
                       m2=mean(y2),sd2=sd(y2),.N),list(j1,x)]
  archive.put(cstatsx=cstats,m="mean and sd by cluster for stayers by x",file = arch_static)

  mini_model = model.mini2.estimate(jdata,sdata,norm = 4,fixb=T,withx=T)
  sdata.sim  = model.mini2.impute.stayers(sdata,mini_model)
  mini_model$vdec = lin.projax(sdata.sim,"y1_imp","k_imp","j1")

  archive.put(mini_model_x=mini_model,m="main mini model estimation with X",file = arch_static)
}


server.static.mini.estimate.main.randomclusters <- function(){
  arch_static = "res-2003-static.dat"

  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  clus = unique(sdata[,list(fid=f1,clus=j1)])
  clus[,clus:=sample.int(10,.N,replace=T)]
  jdata =  cluster.append.data(jdata,clus)
  sdata =  cluster.append.data(sdata,clus)

  mstats = jdata[,list(m1=mean(y1),sd1=sd(y1),
                       m2=mean(y2),sd2=sd(y2),.N),list(j1,j2)]
  cstats = sdata[,list(m1=mean(y1),sd1=sd(y1),
                       m2=mean(y2),sd2=sd(y2),.N),list(j1)]

  mini_model = model.mini2.estimate(jdata,sdata,norm = 4,fixb=T)
  sdata.sim  = model.mini2.impute.stayers(sdata,mini_model)
  vdec_minimodel = lin.proja(sdata.sim,"y1_imp","k_imp","j1")
}



#' estimate the mini-model once, then resimulate from it, then re-cluster
#' and re-estimate. repeat theis many times
server.static.mini.estimate.bootstrap.leave10out <- function() {

  arch = "res-2003-static.dat"
  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))

  fids = sdata[,.N,f1][N>=25,f1]
  sdata = sdata[f1 %in% fids][f2 %in% fids]
  jdata = jdata[f1 %in% fids][f2 %in% fids]

  # reclusterleaving 10% out
  wids.mo = jdata[,wid]
  wids.s1 = sdata[!(wid%in%wids.mo),list(wid=sample(wid,ceiling(.N*0.9))),f1][,wid]
  #wids.s1 = union(
  #  sdata[!(wid%in%wids.mo),list(wid=sample(wid,1)),f1][,wid],
  #  sdata[!(wid%in%wids.mo),list(wid=sample(wid,ceiling(0.9*.N)))][,wid])

  cdata = sdata[(wid%in%wids.s1),list(fid=f1,lw=y1)]
  #cdata = sdata[,list(fid=f1,lw=y1)]
  clus  = cluster.firms.data(cdata,ncluster=10,nw=40,nstart=300,step=50)
  sdata = cluster.append.data(sdata,clus$clus)
  jdata = cluster.append.data(jdata,clus$clus)

  # estimate using 10% left out
  wids.s2 = sdata[(wid %in% wids.s1),sample(wid,ceiling(0.1/0.9*.N))]
  mini_model_true = model.mini2.estimate(jdata,sdata[(wid %in% wids.s2)],norm = 4,fixb=T)
  mini_model_true = model.mini2.estimate(jdata,sdata[!(wid %in% wids.s1)],norm = 4,fixb=T)
  # mini_model_true = model.mini2.estimate(jdata,sdata,norm = 4,fixb=T)

  sdata.sim       = model.mini2.impute.stayers(ssample(sdata,0.1),mini_model_true)
  vdec_minimodel  = lin.proja(sdata.sim,"y1_imp","k_imp","j1")

  rrr=data.frame()
  rr_mix = list()
  rr_mini = list()
  for (i in 1:100) {
    sdata.sim = model.mini2.impute.stayers(sdata,mini_model_true)
    jdata.sim = model.mini2.impute.movers(jdata,mini_model_true)
    sdata.sim[,y1:=y1_imp]
    sdata.sim[,y2:=y2_imp]
    jdata.sim[,y1:=y1_imp]
    jdata.sim[,y2:=y2_imp]
    vdec_dir_100_bc  = lin.proja(ssample(sdata.sim,0.1),"y1","k_imp","j1")

    # take some stayers (not movers) out from clustering
    wids.mo = jdata.sim[,wid]

    # sample overall, to avoid treating small firms differently
    wids.s1 = union(
                sdata.sim[!(wid%in%wids.mo),list(wid=sample(wid,1)),f1][,wid],
                sdata.sim[!(wid%in%wids.mo),list(wid=sample(wid,ceiling(0.9*.N)))][,wid])
    vdec_dir_90_bc  = lin.proja(ssample(sdata.sim[wid%in%wids.s1],0.1),"y1","k_imp","j1")
    vdec_dir_10_bc  = lin.proja(sdata.sim[!(wid%in%wids.s1)],"y1","k_imp","j1")

    # estimate LIML using true clusters
    mini_model_bis   = model.mini2.estimate(jdata.sim,sdata.sim,norm = 4,fixb=T)
    sdata.sim.bis    = model.mini2.impute.stayers(ssample(sdata.sim,0.1),mini_model_bis)
    vdec_liml_100_bc = lin.proja(sdata.sim.bis,"y1_imp","k_imp","j1")

    mini_model_bis  = model.mini2.estimate(jdata.sim,sdata.sim[(wid %in% wids.s1)],norm = 4,fixb=T)
    sdata.sim.bis   = model.mini2.impute.stayers(ssample(sdata.sim,0.1),mini_model_bis)
    vdec_liml_90_bc = lin.proja(sdata.sim.bis,"y1_imp","k_imp","j1")

    mini_model_bis  = model.mini2.estimate(jdata.sim,sdata.sim[!(wid %in% wids.s1)],norm = 4,fixb=T)
    sdata.sim.bis   = model.mini2.impute.stayers(ssample(sdata.sim,0.1),mini_model_bis)
    vdec_liml_10_bc = lin.proja(sdata.sim.bis,"y1_imp","k_imp","j1")

    # --- recluster ---- #
    cdata = sdata.sim[(wid%in%wids.s1),list(fid=f1,lw=y1)]
    clus  = cluster.firms.data(cdata,ncluster=10,nw=40,nstart=300,step=50)
    sdata.sim = cluster.append.data(sdata.sim,clus$clus)
    jdata.sim = cluster.append.data(jdata.sim,clus$clus)
    vdec_dir_90_ac  = lin.proja(ssample(sdata.sim[wid%in%wids.s1],0.1),"y1","k_imp","j1")
    vdec_dir_10_ac  = lin.proja(sdata.sim[!(wid%in%wids.s1)],"y1","k_imp","j1")
    vdec_dir_100_ac  = lin.proja(ssample(sdata.sim,0.1),"y1","k_imp","j1")

    # estimate mini-model
    mini_model_bis = model.mini2.estimate(jdata.sim,sdata.sim[!(wid %in% wids.s1)],norm = 4,fixb=T)
    sdata.sim.bis  = model.mini2.impute.stayers(ssample(sdata.sim,0.1),mini_model_bis)
    vdec_liml_10_ac = lin.proja(sdata.sim.bis,"y1_imp","k_imp","j1")
    rr_mini[[paste(i)]] = mini_model_bis

    mini_model_bis = model.mini2.estimate(jdata.sim,sdata.sim[(wid %in% wids.s1)],norm = 4,fixb=T)
    sdata.sim.bis  = model.mini2.impute.stayers(ssample(sdata.sim,0.1),mini_model_bis)
    vdec_liml_90_ac = lin.proja(sdata.sim.bis,"y1_imp","k_imp","j1")

    mini_model_bis = model.mini2.estimate(jdata.sim,sdata.sim,norm = 4,fixb=T)
    sdata.sim.bis  = model.mini2.impute.stayers(ssample(sdata.sim,0.1),mini_model_bis)
    vdec_liml_100_ac = lin.proja(sdata.sim.bis,"y1_imp","k_imp","j1")


    rr2 = rbind(
      data.frame(vdec_liml_90_ac$stats ,where="vdec_liml_90_ac"),
      data.frame(vdec_liml_10_ac$stats ,where="vdec_liml_10_ac"),
      data.frame(vdec_dir_10_ac$stats  ,where="vdec_dir_10_ac"),
      data.frame(vdec_dir_90_ac$stats  ,where="vdec_dir_90_ac"),
      data.frame(vdec_dir_100_ac$stats ,where="vdec_dir_100_ac"),
      data.frame(vdec_liml_10_bc$stats ,where="vdec_liml_10_bc"),
      data.frame(vdec_liml_90_bc$stats ,where="vdec_liml_90_bc"),
      data.frame(vdec_liml_100_bc$stats,where="vdec_liml_100_bc"),
      data.frame(vdec_dir_10_bc$stats  ,where="vdec_dir_10_bc"),
      data.frame(vdec_dir_90_bc$stats  ,where="vdec_dir_90_bc"),
      data.frame(vdec_dir_100_bc$stats ,where="vdec_dir_100_bc"))
    rr2$model="mini"
    rr2$i= i

    rrr = rbind(rrr,rr2)
    save(rrr,rr_mix,rr_mini,file="L:\\Tibo\\qtrdata\\tmp-static-bootstrap-lo2.dat")
    catf("don with rep %i\n",i)
    print(data.table(rrr)[,lapply(.SD,mean),where,.SDcols = c("cor_kl","cov_kl","var_k","var_l","rsq")])
  }

  rrr = data.table(rrr)
  rbind(
    rrr[,lapply(.SD,mean), ,.SDcols = c("cor_kl","cov_kl","var_k","var_l","rsq")],
    rrr[,lapply(.SD,quantile,0.025), ,.SDcols = c("cor_kl","cov_kl","var_k","var_l","rsq")],
    rrr[,lapply(.SD,quantile,0.975), ,.SDcols = c("cor_kl","cov_kl","var_k","var_l","rsq")])


  archive.put(mini_bs_lo = rr_mini ,m="all realisation of the mini model in each bootstrap result",file = arch)
  archive.put(mini_bs_lo_vardec=rrr,m="all variaance decomposition for each boostrap repetition",file = arch)
}


#' estimate the mini-model once, then resimulate from it, then re-cluster
#' and re-estimate. repeat theis many times
server.static.mini.estimate.bootstrap.firmsize <- function() {

  fsize=1

  arch = "res-2003-static.dat"
  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))

  fids = sdata[,.N,f1][N>=fsize,f1]
  sdata = sdata[f1 %in% fids][f2 %in% fids]
  jdata = jdata[f1 %in% fids][f2 %in% fids]

  cdata = sdata[(wid%in%wids.s1),list(fid=f1,lw=y1)]
  clus  = cluster.firms.data(cdata,ncluster=10,nw=40,nstart=300,step=50)
  sdata = cluster.append.data(sdata,clus$clus)
  jdata = cluster.append.data(jdata,clus$clus)

  # estimate using 10% left out
  mini_model_true = model.mini2.estimate(jdata,sdata,norm = 4,fixb=T)
  sdata.sim       = model.mini2.impute.stayers(ssample(sdata,0.1),mini_model_true)
  vdec_minimodel  = lin.proja(sdata.sim,"y1_imp","k_imp","j1")

  rrr=data.frame()
  rr_mix = list()
  rr_mini = list()
  for (i in 1:100) {
    sdata.sim = model.mini2.impute.stayers(sdata,mini_model_true)
    jdata.sim = model.mini2.impute.movers(jdata,mini_model_true)
    sdata.sim[,y1:=y1_imp]
    sdata.sim[,y2:=y2_imp]
    jdata.sim[,y1:=y1_imp]
    jdata.sim[,y2:=y2_imp]
    vdec_dir_100_bc  = lin.proja(ssample(sdata.sim,0.1),"y1","k_imp","j1")

    # estimate LIML using true clusters
    mini_model_bis   = model.mini2.estimate(jdata.sim,sdata.sim,norm = 4,fixb=T)
    sdata.sim.bis    = model.mini2.impute.stayers(ssample(sdata.sim,0.1),mini_model_bis)
    vdec_liml_100_bc = lin.proja(sdata.sim.bis,"y1_imp","k_imp","j1")

    # --- recluster ---- #
    cdata = sdata.sim[,list(fid=f1,lw=y1)]
    clus  = cluster.firms.data(cdata,ncluster=10,nw=40,nstart=300,step=50)
    sdata.sim = cluster.append.data(sdata.sim,clus$clus)
    jdata.sim = cluster.append.data(jdata.sim,clus$clus)
    vdec_dir_100_ac  = lin.proja(ssample(sdata.sim,0.1),"y1","k_imp","j1")

    # estimate mini-model
    mini_model_bis = model.mini2.estimate(jdata.sim,sdata.sim,norm = 4,fixb=T)
    sdata.sim.bis  = model.mini2.impute.stayers(ssample(sdata.sim,0.1),mini_model_bis)
    vdec_liml_100_ac = lin.proja(sdata.sim.bis,"y1_imp","k_imp","j1")
    rr_mini[[paste(i)]] = mini_model_bis

    rr2 = rbind(
      data.frame(vdec_liml_100_ac$stats ,where="vdec_liml_100_ac"),
      data.frame(vdec_dir_100_ac$stats ,where="vdec_dir_100_ac"),
      data.frame(vdec_liml_100_bc$stats,where="vdec_liml_100_bc"),
      data.frame(vdec_dir_100_bc$stats ,where="vdec_dir_100_bc"))
    rr2$model="mini"
    rr2$i= i

    rrr = rbind(rrr,rr2)
    save(rrr,rr_mix,rr_mini,file=sprintf("L:\\Tibo\\qtrdata\\tmp-static-bootstrap-fsize%i.dat",fsize))
    catf("don with rep %i\n",i)
    print(data.table(rrr)[,lapply(.SD,mean),where,.SDcols = c("cor_kl","cov_kl","var_k","var_l","rsq")])
  }

  load(sprintf("L:\\Tibo\\qtrdata\\tmp-static-bootstrap-fsize%i.dat",5))
  archive.put(mini_bs_sf_5 = rr_mini ,m="minimodel bootstrap fsize 5",file = arch_static)
  archive.put(mini_bs_sf_5_vardec=rrr,m="minimodel bootstrap fsize 5 variance decomposition",file = arch_static)
  load(sprintf("L:\\Tibo\\qtrdata\\tmp-static-bootstrap-fsize%i.dat",10))
  archive.put(mini_bs_sf_10 = rr_mini ,m="minimodel bootstrap fsize 10",file = arch_static)
  archive.put(mini_bs_sf_10_vardec=rrr,m="minimodel bootstrap fsize 10 variance decomposition",file = arch_static)
  load(sprintf("L:\\Tibo\\qtrdata\\tmp-static-bootstrap-fsize%i.dat",50))
  archive.put(mini_bs_sf_50 = rr_mini ,m="minimodel bootstrap fsize 50",file = arch_static)
  archive.put(mini_bs_sf_50_vardec=rrr,m="minimodel bootstrap fsize 50 variance decomposition",file = arch_static)

  rrr = data.table(rrr)
  rbind(
    rrr[,lapply(.SD,mean), ,.SDcols = c("cor_kl","cov_kl","var_k","var_l","rsq")],
    rrr[,lapply(.SD,quantile,0.025), ,.SDcols = c("cor_kl","cov_kl","var_k","var_l","rsq")],
    rrr[,lapply(.SD,quantile,0.975), ,.SDcols = c("cor_kl","cov_kl","var_k","var_l","rsq")])
}



#' estimate the mini-model with different cluster sizes
server.static.mini.estimate.ksize <- function() {

  arch = "res-2003-static.dat"
  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))

  load(file="L:\\Tibo\\qtrdata\\tmp-2003-clusters.dat")
  rr_models = list()
  rr_vardecs =list()
  for (k in 3:20) {
    jdata =  cluster.append.data(jdata,res[[paste(k)]])
    sdata =  cluster.append.data(sdata,res[[paste(k)]])
    mini_model = model.mini2.estimate(jdata,sdata,norm = 4,fixb=T)
    sdata.sim  = model.mini2.impute.stayers(sdata,mini_model)
    vdec_minimodel = lin.proja(sdata.sim,"y1_imp","k_imp","j1")
    vdec_minimodel$between_var_explained = sdata[, list(mean(y1),.N),j1][,wtd.var(V1,N)]/sdata[, list(mean(y1),.N),f1][,wtd.var(V1,N)]*100
    rr_models[[paste(k)]] = mini_model
    rr_vardecs[[paste(k)]] = vdec_minimodel
  }

  archive.put(mini_ksize = rr_models ,m="mini model for different k sizes",file = arch)
  archive.put(mini_ksize_vardec=rr_vardecs,m="all variance decomposition for each cluster size",file = arch)
}


# running on full data
server.static.mini.estimate.fulldata <- function() {
  arch_static = "res-2003-static.dat"

  load("L:\\Tibo\\qtrdata\\tmp-2003-allc.dat")
  cdata = sdata[,list(fid=f1,lw=y1)]
  clus  = cluster.firms.data(cdata,ncluster=10,nw=40,nstart=500,step=100)
  sdata = cluster.append.data(sdata,clus)
  jdata = cluster.append.data(jdata,clus)
  jdata = jdata[!is.na(j2)]
  sdata = sdata[!is.na(j2)]

  mini_model = model.mini2.estimate(jdata,sdata,norm = 4,fixb=T)
  sdata.sim  = model.mini2.impute.stayers(sdata,mini_model)
  vdec_minimodel = lin.proja(sdata.sim,"y1_imp","k_imp","j1")

  archive.put(mini_model_full=mini_model,m="mini model estimated on every one",file = arch_static)
  archive.put(mini_model_full_vardec=vdec_minimodel,m="mini model estimated on every one, variance decomposition",file = arch_static)
}

server.static.mini.estimate.splits <- function() {
  arch_static = "res-2003-static.dat"

  nc = 8 # nuber of clusters
  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,ageg:= (age<=30) + 2*((age >= 31)&(age<=50)) + 3*(age>50)]
  sdata[,x := educ + 3*(ageg-1)]
  sdata[, x := as.integer(x)]

  clusAndEst = function(jdata,sdata,nc=10) {
      cdata = sdata[,list(fid=f1,lw=y1)]
      clus  = cluster.firms.data(cdata,ncluster=nc,nw=40,nstart=300,step=50)
      sdata = cluster.append.data(sdata,clus$clus)
      jdata = cluster.append.data(jdata,clus$clus)

      jdata = jdata[!is.na(j1*j2)]
      sdata = sdata[!is.na(j1*j2)]

      mini_model = model.mini2.estimate(jdata,sdata,norm = 4,fixb=T)
      sdata.sim  = model.mini2.impute.stayers(sdata,mini_model)
      mini_model$vdec = lin.proja(sdata.sim,"y1_imp","k_imp","j1")
      return(mini_model)
  }

  # select one industry
  inds = sdata[,unique(ind1)]
  for (ind in inds) {
    fids = sdata[ind1==ind,unique(f1)]
    mini1 = clusAndEst(jdata[(f1 %in% fids)&(f2 %in% fids)],sdata[(f1 %in% fids)&(f2 %in% fids)],nc=nc)
    res[[sprintf("ind-%s-c%i",ind,nc)]] = mini1
  }

  # split by educ
  educs = sdata[,unique(educ)]
  for (e in educs) {
    wids = sdata[educ==e,unique(wid)]
    mini1 = clusAndEst(jdata[wid %in% wids],sdata[wid %in% wids],nc=nc)
    res[[sprintf("educ-%i-c%i",e,nc)]] = mini1
  }

  # split by age
  sdata[,acat := cut(age,breaks = c(0,30,50,150))]
  acats = sdata[,unique(acat)]
  for (ag in acats) {
      wids = sdata[acat==ag,unique(wid)]
      mini1 = clusAndEst(jdata[wid %in% wids],sdata[wid %in% wids],nc=nc)
      res[[sprintf("age-%s-c%i",ag,nc)]] = mini1
  }

  archive.put(mini_splits=res,m="mini model for different splits based on observables",file = arch_static)
}

server.static.mini.d2003.estimate.within_re <- function() {

  rr = data.frame()

  for (nc in (5:30)) {

    load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
    sdata[,x:=1][,y1_bu:=y1]
    sdata = sdata[move==FALSE]

    # on connected set
    #f0s   = get.largest.conset.fid(sim$jdata)
    #jdata = sim$jdata[f1%in%f0s][f2%in%f0s]
    #sdata = sim$sdata[f1%in%f0s]

    sim = list(jdata=jdata,sdata=sdata)
    rm(sdata,jdata)

    ms    = grouping.getMeasures(sim,"ecdf",Nw=20,y_var = "y1")
    grps  = grouping.classify.once(ms,k = nc,nstart = 1000,iter.max = 200,step=100)
    sim   = grouping.append(sim,grps$best_cluster)

    # --- MINI ESTIMATION ---- #
    model_minilin = m2.mini.estimate(sim$jdata,sim$sdata,method="linear")
    # --- WITHIN RE ESTIMATION ---- #
    sim$jdata[,psi1 := model_minilin$A1[j1]]
    sim$jdata[,psi2 := model_minilin$A2[j2]]
    S = sim$sdata[,.N,j1][order(j1),N]

    rrr = list()
    rrr$nc = nc
    rrr$y_var       = model_minilin$vdec$cc[1,1]
    rrr$cluster_var = model_minilin$vdec$cc[3,3]

    omega = m2.trace.reExt(sim$jdata)
    rrr$omega_var = wt.mean(omega,S)
    rrr$omega_var_pos = wt.mean(pmax(omega,0),S)

    omega = m2.trace.reExt(sim$jdata,2)
    rrr$omega_var_sub2 = wt.mean(omega,S)
    rrr$omega_var_sub2_pos = wt.mean(pmax(omega,0),S)

    rr = rbind(rr,data.frame(rrr))
    print(rr)
  }

  # 20  clusters  4.07% ( 2.1% within )
  # 20 (sub=2)    4.54% ( 2.6% within )
  # 10 clusters   4.69%
  # 5  clusters   3.50% ( 2.1% within )
  # 5 (sub=2)     4.2%  ( 2.7% within )

}

# ==== MIXTURE MODEL =====

server.static.mixture.d2003.estimate <- function() {

  sim = server.static.data()
  grps  = res.load("m2-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  # --- MAIN ESTIMATION ---- #
  set.seed(87954352)

  # --- use cluster ---#
  cl = makeCluster(local_opts$number_of_clusters)
  clusterEvalQ(cl,require(blmrep))

  # --- MAIN ESTIMATION -  STATIONARY INTERACTIONS ---- #
  ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=50,est_nbest=10,sdata_subsample=0.1)
  res_mixt = m2.mixt.estimate.all(sim,nk=6,ctrl,cl)
  res.save("m2-mixt-d2003-main-fixb",res_mixt)

  # --- MAIN ESTIMATION - NON STATIONARY ---- #
  ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=FALSE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=50,est_nbest=10,sdata_subsample=0.1)
  res_mixt = m2.mixt.estimate.all(sim,nk=6,ctrl,cl)
  res.save("m2-mixt-d2003-main-ns",res_mixt)

  stopCluster(cl)
}

server.static.mixture.d2003.estimate.withx <- function() {

  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,ageg:= (age<=30) + 2*((age >= 31)&(age<=50)) + 3*(age>50)]
  sdata[,x := educ + 3*(ageg-1)]
  sdata[,x := as.integer(x)]
  sdata[,y1_bu:=y1]
  # we include the m
  sim = list(jdata=jdata,sdata=sdata)
  rm(sdata,jdata)

  # get clsuters
  grps  = rkiv0.load("m2-mixt-y2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  # get main estimate
  res_main = rkiv0.load("m2-mixt-y2003-main-fixb")

  # estimate using observables
  ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=50,est_nbest=10,sdata_subsample=0.1)

  sim$sdata[,sample := rank(runif(.N))/.N<=ctrl$sdata_subsample,j1]
  res_main$model$pk0 = spread(rdim(res_main$model$pk0,10,6),1,9)
  res_main$model = m2.mixt.stayers(sim$sdata[sample==1],res_main$model,ctrl = em.control(ctrl,textapp="stayers"))

  sprop = sim$sdata[,.N,list(x,j1)]

  rs    = rkiv0.start("m2-mixt-d2003-main-withx",
                      info="static 2003 main estimate, with observables")
  rkiv0.put(rs,list(mix=res_main,sprop=sprop))

}


# here we retain 10% of the stayers in classification
server.static.mixture.d2003.estimate.holdout <- function() {

  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]
  sdata = sdata[move==FALSE]
  sim = list(jdata=jdata,sdata=sdata)
  rm(sdata,jdata)

  rs    = rkiv0.start("m2-mixt-d2003-main-holdout",
                      info="static 2003 main estimate, holding stayers in clustering")

  # create holdout
  sim$sdata[,HO := rank(runif(.N))/.N <= 0.1 ,f1]

  ms    = grouping.getMeasures(list(sdata=sim$sdata[HO==FALSE],jdata=sim$jdata),"ecdf",Nw=20,y_var = "y1")
  grps  = grouping.classify.once(ms,k = 10,nstart = 1000,iter.max = 200,step=100)
  sim   = grouping.append(sim,grps$best_cluster)

  ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=34,est_nbest=10,sdata_subsample=1,sdata_subredraw=FALSE)
  cl = makeCluster(17)
  clusterEvalQ(cl,require(blm))

  sim$sdata[,sample:=HO]
  res_mixt = m2.mixt.estimate.all(sim,nk=6,ctrl,cl)

  stopCluster(cl)
  rkiv0.put(rs,res_mixt)

}

server.static.mixture.d2003.estimate.holdout2 <- function() {

  # try to include half the stayers in clustering
  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]

  # split movers in 2
  jdata[,splitdata := rank(runif(.N)) - 0.5 +0.1*(runif(1)-0.5) <= .N/2,f1]
  mwids = jdata[splitdata==TRUE,unique(wid)]
  sdata = sdata[!(wid %in% mwids)]
  jdata = jdata[wid %in% mwids]
  sim = list(jdata=jdata,sdata=sdata)
  rm(sdata,jdata)

  rs    = rkiv0.start("m2-mixt-d2003-main-holdout2",
                      info="static 2003 main estimate, holding stayers in clustering")

  ms    = grouping.getMeasures(sim,"ecdf",Nw=20,y_var = "y1")
  grps  = grouping.classify.once(ms,k = 10,nstart = 1000,iter.max = 200,step=100)
  sim   = grouping.append(sim,grps$best_cluster)

  ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=50,est_nbest=10,sdata_subsample=0.1,sdata_subredraw=TRUE)
  cl = makeCluster(17)
  clusterEvalQ(cl,require(blm))
  res_mixt = m2.mixt.estimate.all(list(sdata=sim$sdata[move==FALSE],jdata=sim$jdata),nk=6,ctrl,cl)

  stopCluster(cl)
  rkiv0.put(rs,res_mixt)
}

server.static.mixture.d2003.fit <- function() {

  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]
  sim = list(jdata=jdata,sdata=sdata)
  rm(sdata,jdata)

  # get clsuters
  grps  = rkiv0.load("m2-mixt-y2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  # get main estimate
  res_main = rkiv0.load("m2-mixt-y2003-main-fixb")

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

  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,x:=1]
  sdata=sdata[move==FALSE]
  sim = list(jdata=jdata,sdata=sdata)
  rm(sdata,jdata)

  set.seed(87954352)
  rs    = rkiv0.start("m2-mixt-y2003-changeK",
                      info="static 2003 estimation with different number of clusters")
  ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=45,est_nbest=10,sdata_subsample=0.1)

  cl = makeCluster(15)
  clusterEvalQ(cl,require(blm))

  rr_mixt = list()
  for (nf_size in c(3:15,17,20)) {
    tryCatch({
      ms    = grouping.getMeasures(sim,"ecdf",Nw=20,y_var = "y1")
      grps  = grouping.classify.once(ms,k = nf_size,nstart = 500,iter.max = 200)
      sim   = grouping.append(sim,grps$best_cluster)
      res_mixt = m2.mixt.estimate.all(sim,nk=6,ctrl,cl=cl)
      rr_mixt[[paste("nf",nf_size,sep="-")]] = res_mixt
    })}
  stopCluster(cl)
  rkiv0.put(rs,rr_mixt)
}

#' this estimates the model with different values of K to check
#' robustness. We recluster every time
server.static.mixture.estimate.robust.nk <- function() {

  # load data
  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,x:=1]
  sim = list(jdata=jdata,sdata=sdata)
  rm(sdata,jdata)

  # get clsuters
  grps  = rkiv0.load("m2-mixt-y2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  rs    = rkiv0.start("m2-mixt-d2003-change-nk",
                      info="static 2003 estimation with different number of worker types")
  ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=45,est_nbest=5,sdata_subsample=0.1)

  cl = makeCluster(15)
  clusterEvalQ(cl,require(blm))

  rr_mixt = list()
  for (nk_size in c(3:5,7:9)) {
    tryCatch({
      res_mixt = m2.mixt.estimate.all(sim,nk=nk_size,ctrl=ctrl,cl=cl)
      res_mixt$second_stage_reps_all=NA
      rr_mixt[[paste("nk",nk_size,sep="-")]] = res_mixt
    })}
  stopClsuter(cl)
  rkiv0.put(rs,rr_mixt)

  rr = data.table(ldply(rr_mixt,function(r) {
    dd = data.frame(r$vdec$stats)
    dd$nk = r$model$nk
    dd$con = r$connectedness
    dd
  }))

  rr[order(nk)]
}


#' we estimate for small firms and for large firms
server.static.mixture.estimate.robust.fsize <- function() {

  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,x:=1]
  sdata=sdata[move==FALSE]
  sim = list(jdata=jdata,sdata=sdata)
  rm(sdata,jdata)

  # select firms larger or smaller than 50
  fids = unique(sim$sdata[,.N,f1][N>50, f1])

  sim$sdata = sim$sdata[f1 %in% fids]
  sim$jdata = sim$jdata[f1 %in% fids][f2 %in% fids]

  set.seed(87954352)
  rs    = rkiv0.start("m2-mixt-d2003-firmsize",
                      info="estimating conditioning on different firm sizes")
  ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=45,est_nbest=10,sdata_subsample=0.1)

  cl = makeCluster(15)
  clusterEvalQ(cl,require(blm))
  ms    = grouping.getMeasures(sim,"ecdf",Nw=20,y_var = "y1")
  grps  = grouping.classify.once(ms,k = 10,nstart = 1000,iter.max = 200)
  sim   = grouping.append(sim,grps$best_cluster)
  res_mixt = m2.mixt.estimate.all(sim,nk=6,ctrl,cl=cl)
  stopCluster(cl)

  save(res_mixt,file="tmp-size50.dat")
  res_mixt_bt_50 = res_mixt
  load("tmp-size50.dat")

  rkiv0.put(rs,list(mixt_leq_50=res_mixt,mixt_g_50=res_mixt_bt_50))
}



#' Boostrapping the main results. We use the value for the bias-corrected
#' estimates.
server.static.mixture.estimate.boostrap <- function(){

  sim = server.static.data()
  grps  = res.load("m2-mixt-d2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  # get main estimate
  res_main  = res.load("m2-mixt-y2003-main-fixb")
  model_true = res_main$model

  ctrl      = em.control(nplot=1000,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=45,est_nbest=10,sdata_subsample=0.1)

  cl = makeCluster(15)
  clusterEvalQ(cl,require(blm))

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
  rs      = rkiv0.start("m2-mixt-d2003-bootstrap-leg2")
  res.save("m2-mixt-d2003-bootstrap",list(vdecs=rrr,mixt_all=rr_mixt2))
}


#' Boostrapping the fixed-effect estimator
server.static.mixture.estimate.boostrap.akm <- function(){

  # get data
  arch_static = "res-2003-static.dat"
  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,x:=1]
  sdata = sdata[move==FALSE]
  sim = list(sdata=sdata,jdata=jdata)

  # get groups
  grps  = rkiv0.load("m2-mixt-y2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  # get main estimate
  res_main  = rkiv0.load("m2-mixt-y2003-main-fixb")
  model_true = res_main$model

  # now we bootstrap with clustering
  rrr=data.frame()
  rr_mixt = list()
  for (i in reps) {
    tryCatch({
      # we use the regression coef estimated in the dynamic
      sdata.sim = m2.mixt.impute.stayers(sim$sdata,model_true,rho21s=0.681,stationary=TRUE)
      jdata.sim = m2.mixt.impute.movers(sim$jdata,model_true)
      sdata.sim[,y1:=y1_imp][,y2:=y2_imp][,jt:=j1][,y1_bu:=y1]
      jdata.sim[,y1:=y1_imp][,y2:=y2_imp]
      sim.sp = list(jdata=jdata.sim,sdata=sdata.sim)

      # --- estimate AKM --- #
      akm_cov = m2.fe.all(sim.sp)

      rr2 = data.frame(minib$vdec$stats)

      rr2$rep=i
      rrr = rbind(rrr,rr2)
      flog.info("done with rep %i",i)

      print(rrr)
    }, error = function(e) {catf("error in boot strap rep %i!\n",i);print(e)})
  }

  stopCluster(cl)

  load("L:\\Tibo\\qtrdata\\tmp-static-bootstrap-mixt-45s-100to200.dat")
  rr_mixt2 = lapply(rr_mixt,function(r) { r$second_stage_reps_all=NULL;r})
  rs      = rkiv0.start("m2-mixt-d2003-bootstrap-leg2")
  rs$info = "parametric bootstrap of static mixture on 2003 data reps 100 to 200"
  rkiv0.put(rs,list(vdecs=rrr,mixt_all=rr_mixt2))



  load("L:\\Tibo\\qtrdata\\tmp-static-bootstrap-mixt.dat")
  archive.put(mixt_bs = rr_mixt ,m="all realisation of the mixture model in each bootstrap result",file = arch_static)
}

server.static.estimate.poachingrank <- function() {

  arch_static = "res-2003-static.dat"
  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)
  #sim$jdata$fids=NULL

  frank <- sim$sdata[,mean(y1),f1][,list(fq=ceiling(3*rank(V1)/.N),f1)]
  setkey(frank,f1)
  setkey(sim$jdata,f1)
  sim$jdata[,fq1 := frank[sim$jdata,fq]]
  setkey(sim$jdata,f2)
  sim$jdata[,fq2 := frank[sim$jdata,fq]]

  # create the measure
  # for each firm we count number of movers coming/going going to each firm-quartile
  rr1=sim$jdata[,list(.N,dir="out"),list(f=f1,fq=fq2)]
  rr2=sim$jdata[,list(.N,dir="in"),list(f=f2,fq=fq1)]
  rr=rbind(rr1,rr2)
  rr[,S:=sum(N),f]
  rr[,r:=N/S]
  rr2 = rr[,list(f1=f,measure=paste(dir,fq,sep="-"),value=r)]

  # select firms with at least nm movers
  # fids = rr[S>=10,unique(f)]
  # rr = rr[f %in% fids]
  #setkey(rr,f)

  # append additional measure
  rr_more = sim$sdata[f1 %in% unique(rr2$f1),list(wm=mean(y1),wsd=sd(y1)),f1]
  rr_more[is.na(wsd),wsd:=0]
  rr_more = data.table(melt(rr_more,id.vars = "f1") )
  setnames(rr_more,"variable","measure")
  rr = rbind(rr_more,rr2)

  # get size from cross-section
  ddf = sim$sdata[,list(nw=.N),f1]
  setkey(ddf,f1)
  setkey(rr,f1)
  rr = ddf[rr]

  # standardize weighted
  rr[,m0  := sum(nw*value)/sum(nw), list(measure)]
  rr[,sd0 := sqrt(sum(nw*(value-m0)^2/sum(nw))), list(measure)]
  rr[,rf  := (value-m0)/sd0, list(measure)]

  setkey(rr,f1)
  M = acast(rr[str_detect(measure,"wm|wsd")],f1 ~ measure,fill=0,value.var = "rf")
  W = rr[,nw[1],f1][,V1] # use total size

  # compute the groups
  grps  = grouping.classify.once(list(M=M,W=W,Nw=ncol(M),measure="move",N=nrow(M),tsize=0,discrepency=0),k = 10,nstart = 1000,iter.max = 200)
  sim   = grouping.append(sim,grps$best_cluster,drop=T)

  # moving patterns
  acast(sim$jdata[,.N,list(j1,j2)],j1~j2,fill=0)
  sim$sdata[,mean(y1),j1][order(j1)]

  # between firm variance
  my = sim$sdata[,list(m=mean(y1),.N),f1][,sum(N*m)/sum(N)]
  v1 = sim$sdata[,list(m=mean(y1),.N),f1][,sum(N*(m-my)^2)/sum(N)]
  my = sim$sdata[,list(m=mean(y1),.N),j1][,sum(N*m)/sum(N)]
  v2 = sim$sdata[,list(m=mean(y1),.N),j1][,sum(N*(m-my)^2)/sum(N)]
  v2/v1

  my = sim$jdata[,list(m=mean(y1),.N),f1][,sum(N*m)/sum(N)]
  v1 = sim$jdata[,list(m=mean(y1),.N),f1][,sum(N*(m-my)^2)/sum(N)]
  my = sim$jdata[,list(m=mean(y1),.N),j1][,sum(N*m)/sum(N)]
  v2 = sim$jdata[,list(m=mean(y1),.N),j1][,sum(N*(m-my)^2)/sum(N)]
  v2/v1
}


server.static.estimate.clustersplits <- function() {
  # === EXERCICE 2 ===
  # we take our clusters, we then split each cluster into
  # 2 sub-groups, one with high IN, one with low IN
  # then use it for estimation
  arch_static = "res-2003-static.dat"
  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)

  frank <- sim$sdata[,mean(y1),f1][,list(fq=rank(V1)/.N,f1)]
  setkey(frank,f1)
  setkey(sim$jdata,f1)
  sim$jdata[,fq1 := frank[sim$jdata,fq]]
  setkey(sim$jdata,f2)
  sim$jdata[,fq2 := frank[sim$jdata,fq]]

  # create the measure
  # for each firm we count number of movers coming/going going to each firm-quartile
  if (measure=="poachingrk") {
    rr = sim$jdata[,list(r=median(fq1)),list(f=f2)]

  } else if (measure=="poachingrk-nmover") {
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
  grps  = rkiv0.load("m2-mixt-y2003-groups")
  clus  = data.table(f = names(grps$best_cluster), j = grps$best_cluster)
  setkey(rr,f)
  setkey(clus,f)
  rr = clus[rr]

  # split in 2 within
  rr[,g:=r<median(r),j]
  rr[,jn := j + 10*g]
  clus = rr[,jn]
  names(clus) = rr[,f]

  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)
  sim   = grouping.append(sim,clus,drop=T)

  # we then estimate the model
  ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=45,est_nbest=10,sdata_subsample=0.1,sdata_subredraw=TRUE)
  cl = makeCluster(15)
  clusterEvalQ(cl,require(blm))
  res_mixt = m2.mixt.estimate.all(list(sdata=sim$sdata[move==FALSE],jdata=sim$jdata),nk=6,ctrl,cl)
  stopCluster(cl)

  rs      = rkiv0.start("m2-mixt-d2003-poachingrk-va")
  rs$info = "splitting classes with rank in value added"
  rkiv0.put(rs,res_mixt)
}

server.static.mixt.estimate.reclassify <- function() {

  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]
  sdata= sdata[move==FALSE]
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)

  # get groups
  grps  = rkiv0.load("m2-mixt-y2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  sim$sdata[,jm1:=j1]
  sim$jdata[,jm1:=j1]
  sim$jdata[,jm2:=j2]

  # get main estimate
  res_main  = rkiv0.load("m2-mixt-y2003-main-fixb")
  model = res_main$model

  sim$jdata[,splitdata := rank(runif(.N)) - 0.5 +0.1*(runif(1)-0.5) <= .N/2,f1]

  # we then estimate the model
  ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=45,est_nbest=10,sdata_subsample=0.1,sdata_subredraw=TRUE)
  cl = makeCluster(15)
  clusterEvalQ(cl,require(blm))

  rr_mixt = list()
  for (iter in 1:5) {

    # compute the posterior probability on stayers
    sdata.pos = sim$sdata[,{
      likm = rep(0,model$nf)
      # iterate on firm types
      for (ii in 1:.N) {
        ltau     = log(model$pk0[1,,])
        lnorm1   = lognormpdf(y1[ii], model$A1, model$S1)
        lall     = ltau + lnorm1
        likm     = likm + blm:::logRowSumExp(lall)
      }
      # draw the type of the firm
      list(jp=1:model$nf,lik=likm,N=rep(.N,model$nf),jo=0)
    },list(fid=f1,jt=j1,jm=jm1)]

    # compute the posterior probability on movers
    jdata.pos.p1 = sim$jdata[splitdata==TRUE,{
      likm = rep(0,model$nf)
      # iterate on firm types
      for (ii in 1:.N) {
        ltau     = log(rdim(model$pk1,model$nf,model$nf,model$nk)[,jo,])
        lnorm1   = lognormpdf(y1[ii], model$A1, model$S1)
        lnorm2   = lognormpdf(y2[ii], model$A2[jo,], model$S2[jo,])
        lall     = ltau + lnorm1  + spread(lnorm2,1,model$nf)
        likm     = likm + blm:::logRowSumExp(lall)
      }
      # draw the type of the firm
      list(jp=1:model$nf,lik=likm,N=rep(.N,model$nf))
    },list(fid=f1,jt=j1,jo=j2,jm=jm1)]

    jdata.pos.p2 = sim$jdata[splitdata==TRUE,{
      likm = rep(0,model$nf)
      # iterate on firm types
      for (ii in 1:.N) {
        ltau     = log(rdim(model$pk1,model$nf,model$nf,model$nk)[jo,,])
        lnorm2   = lognormpdf(y2[ii], model$A2, model$S2)
        lnorm1   = lognormpdf(y1[ii], model$A1[jo,], model$S1[jo,])
        lall     = ltau + lnorm2 + spread(lnorm1,1,model$nf)
        likm     = likm + blm:::logRowSumExp(lall)
      }
      # draw the type of the firm
      list(jp=1:model$nf,lik=likm,N=rep(.N,model$nf))
    },list(fid=f2,jt=j2,jo=j1,jm=jm2)]

    # combine all likelihoods
    jdata.pos.p2[,lambda := 1]
    jdata.pos.p1[,lambda := 1]
    sdata.pos[,lambda := 1]
    data.pos  = rbind(sdata.pos,jdata.pos.p1,jdata.pos.p2)
    data.pos  = data.pos[,list(lik=sum(lambda*lik),N=sum(N),Ns=sum(N[jo==0])),list(fid,jp,jt,jm)]
    data.pos  = data.pos[,{i = which.max(lik);list(jp=jp[i],N=N[i],Ns=Ns[i])},list(fid,jt,jm)]

    # correlation
    flog.info("cor=%f cor2=%f ", data.pos[,wt.cor(jt,jp,N)],data.pos[,wt.cor(jp,jm,N)])

    # create clusters out of it
    clus        = data.pos[,jp]
    names(clus) = data.pos[,fid]

    # attach groups
    sim   = grouping.append(sim,clus)

    # we then estimate the model
    ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                           sd_floor=1e-7,posterior_reg=1e-8,
                           est_rep=45,est_nbest=10,sdata_subsample=0.1,sdata_subredraw=TRUE)
    res_mixt = m2.mixt.estimate.all(list(sdata=sim$sdata[move==FALSE],jdata=sim$jdata[splitdata==FALSE]),nk=6,ctrl,cl)
    model    = res_mixt$model
    rr_mixt[[paste(iter)]] = res_mixt
  }

  stopCluster(cl)

  rs    = rkiv0.start("m2-mixt-d2003-main-reclassify-holdout",
                      info="static 2003 estimate, use half of the movers to reclassify firms")
  rkiv0.put(rs,rr_mixt)
}

#' We are going to reclassify using indepenent updating in
#' with the conditional model. Here we will start from the
#' FE classification
server.static.mixt.estimate.dirty_iteration <- function() {
  require(Ckmeans.1d.dp)

  # ===== POACHING RANK ========= #
  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]
  sdata= sdata[move==FALSE]
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)

  load("L:\\Tibo\\qtrdata\\tmp-2003-prank.dat")
  fmean = data_prank[from_j2j>0][from_u>0][,list(from_u/(from_j2j + from_u),N=(from_j2j + from_u)),list(f1 = fid)]
  clusters   = Ckmeans.1d.dp(fmean$V1, 10,fmean$N)
  grp_akm    = clusters$cluster;names(grp_akm) = fmean$f1
  sim        = grouping.append(sim,grp_akm,drop=T)

  res =  m2.mixt.estimate.reclassify(sim,10)
  rs  = rkiv0.start("m2-mixt-d2003-main-reclassify-startwithprank",
                      info="static 2003 dirty iteration, no split")
  rkiv0.put(rs,res)

  # ===== RATIO OF MOVERS ========= #
  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]

  fmean = sdata[,list(V1=sum(move==TRUE)/.N,.N),f1]

  sdata= sdata[move==FALSE]
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)
  clusters   = Ckmeans.1d.dp(fmean$V1, 10,fmean$N)
  grp_akm    = clusters$cluster;names(grp_akm) = fmean$f1
  sim        = grouping.append(sim,grp_akm,drop=T)

  res =  m2.mixt.estimate.reclassify(sim,10)
  rs  = rkiv0.start("m2-mixt-d2003-main-reclassify-startwithratiomover",
                    info="static 2003 dirty iteration, no split")
  rkiv0.put(rs,res)



  # ===== Using residuals ========= #
  # see section on residuals

  # ===== START WITH MAIN RESULTS ========= #
  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]
  sdata= sdata[move==FALSE]
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)

  grps  = rkiv0.load("m2-mixt-y2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  res =  m2.mixt.estimate.reclassify(sim,10)
  rs    = rkiv0.start("m2-mixt-d2003-main-reclassify-mainresults",
                      info="static 2003 dirty iteration, no split")
  rkiv0.put(rs,res)

  # ===== START WITH VALUE ADDED ========= #
  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]
  sdata= sdata[move==FALSE]
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)

  fmean = sim$sdata[,log(va1[1]/size1[1]),f1][is.finite(V1)]
  clusters   = Ckmeans.1d.dp(fmean$V1, 10)
  grp_akm        = clusters$cluster
  names(grp_akm) = fmean$f1
  sim        = grouping.append(sim,grp_akm,drop=T)

  res =  m2.mixt.estimate.reclassify(sim,10)
  rs    = rkiv0.start("m2-mixt-d2003-main-reclassify-startwithva",
                      info="static 2003 dirty iteration, no split")
  rkiv0.put(rs,res)

  # ===== START WITH MEANS ========= #
  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]
  sdata= sdata[move==FALSE]
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)

  fmean = sim$sdata[,list(mean(y1),.N),f1]
  clusters   = Ckmeans.1d.dp(fmean$V1, 10,fmean$N)
  grp_akm    = clusters$cluster; names(grp_akm) = fmean$f1
  sim        = grouping.append(sim,grp_akm,drop=T)

  res =  m2.mixt.estimate.reclassify(sim,20)
  rs    = rkiv0.start("m2-mixt-d2003-main-reclassify-startwithmean",
                      info="static 2003 dirty iteration, no split")
  rkiv0.put(rs,res)

  # ===== START WITH AKM ========= #
  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]
  sdata= sdata[move==FALSE]
  sim=list(sdata=sdata,jdata=jdata)
  rm(sdata,jdata)

  # use AKM grouping
  firm.fe = m2.fe.firms(sim$jdata)

  # classify according to AKM
  require(Ckmeans.1d.dp)
  clusters   = Ckmeans.1d.dp(firm.fe$fe$psi, 10)
  grp_akm        = clusters$cluster
  names(grp_akm) = firm.fe$fe$f1
  sim        = grouping.append(sim,grp_akm,drop=T)

  res = m2.mixt.estimate.reclassify(sim,20)
  rs    = rkiv0.start("m2-mixt-d2003-main-reclassify-startwithakm",
                      info="static 2003 AKM estimate, no split sample")
  rkiv0.put(rs,res)

  ldply(res,function(v) data.frame(v$vdec$stats))


  # combining results
  rr_all = data.frame()
  ll = list(static="m2-mixt-d2003-main-reclassify-mainresults",
            akm="m2-mixt-d2003-main-reclassify-startwithakm",
            mean="m2-mixt-d2003-main-reclassify-startwithmean",
            va = "m2-mixt-d2003-main-reclassify-startwithva",
            prank = "m2-mixt-d2003-main-reclassify-startwithprank",
            mratio = "m2-mixt-d2003-main-reclassify-startwithratiomover",
            resid = "m2-mixt-d2003-main-reclassify-residuals",
            dynanmic = "m4-mixt-d2003-dirty-reclassification")
  for (nn in names(ll)) {
    res = rkiv0.load(ll[[nn]])
    rr = ldply(res,function(v) data.frame(v$vdec$stats))
    rr$start= nn
    rr_all = rbind(rr_all,rr)
  }
  rs    = rkiv0.start("m2-mixt-d2003-main-reclassify-allresults",
                      info="combined results for classification")
  rkiv0.put(rs,rr_all)
}

server.static.mixt.estimate.fit.bs <- function() {

  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]
  sim = list(jdata=jdata,sdata=sdata)
  rm(sdata,jdata)

  # get clsuters
  grps  = rkiv0.load("m2-mixt-y2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  # get main estimate
  res_main = rkiv0.load("m2-mixt-y2003-main-fixb")

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

  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]

  # create some variables
  # get industry for all firms
  firm_info = unique(sdata[,list(f1,ind1)])
  setkey(firm_info,f1)
  setkey(sdata,f2)
  sdata[,ind2:= firm_info[sdata,ind1]]
  sdata[,educ_f := factor(educ)]
  jdata[,educ_f := factor(educ)]

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

  rs    = rkiv0.start("m2-mixt-d2003-residuals-groups",info="static 2003 data cluster outcome for mixture on residuals",location = "serverdata")
  rkiv0.put(rs,grps)
  grps = rkiv0.load("m2-mixt-d2003-residuals-groups")

  # estimation
  ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=50,est_nbest=10,sdata_subsample=0.1)

  cl = makeCluster(17)
  clusterEvalQ(cl,require(blm))
  res_mixt = m2.mixt.estimate.all(sim,nk=6,ctrl,cl)
  stopCluster(cl)

  rs = rkiv0.start("m2-mixt-d2003-main-residuals",info="estimates the mixture model using residual wages")
  rkiv0.put(rs,res_mixt)

  # also do re-classification
  res =  m2.mixt.estimate.reclassify(sim,10)
  rs  = rkiv0.start("m2-mixt-d2003-main-reclassify-residuals",
                    info="static 2003 dirty iteration, no split")
  rkiv0.put(rs,res)



}



# ========== MIXTURE OF MIXTURE MODEL =====

# estimates the mixture of mixture model
server.static.mixture.mixtofmixt <- function() {

  # load data
  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,x:=1]
  sdata <- sdata[move==0]
  sim = list(jdata=jdata,sdata=sdata)
  rm(sdata,jdata)

  # get clsuters
  grps  = rkiv0.load("m2-mixt-y2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  # get main estimate
  res_main = rkiv0.load("m2-mixt-y2003-main-fixb")
  model_mixt = res_main$model
  model_mixt$S1 = 0.9*model_mixt$S1 + 0.1*mean(model_mixt$S1)
  model_mixt$S2 = 0.9*model_mixt$S2 + 0.1*mean(model_mixt$S2)
  model_np   = m2.mixt.np.new.from.ns(model_mixt,3)

  ctrl   = em.control(tol=1e-6,fixb=F,dprior=1.001,maxiter=1000,ncat=1,posterior_reg=1e-7,sd_floor=1e-6)

  # estimate movers
  res_np = m2.mixt.np.movers.estimate(sim$jdata,model_np,ctrl)
  # estimate stayers
  sim$sdata[,sample := rank(runif(.N))/.N<=ctrl$sdata_subsample,j1]
  res_np = m2.mixt.np.stayers.estimate(sim$sdata[sample==1],res_np$model,ctrl)

  rr = blm:::m2.mixt.np.movers.residuals(res_np$model,TRUE,100)

  stayer_share = sim$sdata[,.N]/(sim$sdata[,.N]+sim$jdata[,.N])
  vdec = m2.mixt.np.vdec(res_np$model,nsim = 1e6,stayer_share = stayer_share)
  res_np$vdec = vdec

  rs    = rkiv0.start("m2-mixt-d2003-main-mixtofmixt",
                      info="static 2003 mixture of mixture estimates")
  rkiv0.put(rs,res_np)



}


# ========== PROBABILISTIC APPROACH =====

# this containss function associated with the probabilistic
# evaluation of different models.

server.static.proba.evallik <- function() {

  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]
  sim = list(jdata=jdata,sdata=sdata)
  rm(sdata,jdata)

  # get clsuters
  grps  = rkiv0.load("m2-mixt-y2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)

  # get main estimate
  res_main = rkiv0.load("m2-mixt-y2003-main-fixb")


  rr = m2.proba.importancelik(sim,res_main$model,eps_mix = 0.1,class_draws = 200)
  rr2 = m2.proba.importancelik(sim,res_main$model,eps_mix = 0.5,class_draws = 200)


  rrd = data.table(rr$all[is.finite(rr$lik),])
  rrd[, V := lik - prop_pr + class_prior]
  V = rrd[,V]

  res = data.frame()
  for (l in 1:length(V)) {
    R = rep(0,100)
    for (i in 1:200) {
      vs = sample(V,l,replace = T)
      R[i] = logsumexp(vs) - log(l)
    }
    res = rbind(res,data.frame(m=mean(R),sd=sd(R),l=l))
  }

  res1$eps=0.1
  res2$eps=0.5
  res = rbind(res1,res2)

  ggplot(res,aes(x=l,y=m)) + geom_line() + theme_bw() +
    geom_line(aes(y=m+sd),linetype=2) +
    geom_line(aes(y=m-sd),linetype=2) + scale_x_log10()

  plot(res$m)
}

server.static.proba.feclass <- function() {

  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]
  sdata=sdata[move==FALSE]
  sim = list(jdata=jdata,sdata=sdata)
  rm(sdata,jdata)

  # ru AKM on jdata to extract firm FE
  firm.fe = m2.fe.firms(sim$jdata)

  # classify according to AKM
  require(Ckmeans.1d.dp)
  clusters   = Ckmeans.1d.dp(firm.fe$fe$psi, 10)
  grp_akm        = clusters$cluster
  names(grp_akm) = firm.fe$fe$f1
  sim        = grouping.append(sim,grp_akm,drop=T)
  acast(sim$jdata[,.N,list(j1,j2)],j1~j2,fill=0)

  dd = data.table(f1=names(grp_akm),j=grp_akm)
  dd = merge(firm.fe$fe,dd,by="f1")
  ggplot(dd,aes(x=j,y=psi)) + geom_point() + theme_bw()

  # Estimate BLM model using this classification
  ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=50,est_nbest=10,sdata_subsample=0.1)

  cl = makeCluster(17)
  clusterEvalQ(cl,require(blm))

  res_mixt_akm = m2.mixt.estimate.all(sim,nk=6,ctrl,cl)
  rs    = rkiv0.start("m2-akmgrp-y2003-fixb",info="static 2003 mixture model on group from akm")
  rkiv0.put(rs,res_mixt_akm)

  rr_akm = m2.proba.importancelik(sim,res_mixt_akm$model,eps_mix = 0.1,class_draws = 200)

  # using the same smaple, cluster and estimate
  ms    = grouping.getMeasures(sim,"ecdf",Nw=20,y_var = "y1")
  grps_kmean  = grouping.classify.once(ms,k = 10,nstart = 1000,iter.max = 200,step=100)
  sim   = grouping.append(sim,grps_kmean$best_cluster)
  acast(sim$jdata[,.N,list(j1,j2)],j1~j2,fill=0)

  res_mixt_blm = m2.mixt.estimate.all(sim,nk=6,ctrl,cl)
  rr_blm = m2.proba.importancelik(sim,res_mixt_blm$model,eps_mix = 0.1,class_draws = 200)

  res_akm = prepare(rr_akm)
  res_blm = prepare(rr_blm)
  res_akm$model="akm"
  res_blm$model="blm"
  res = rbind(res_akm,res_blm)

  ggplot(res,aes(x=l,y=m,color=factor(model))) + geom_line() + theme_bw() + scale_x_log10()

  rs    = rkiv0.start("m2-akmgrp-y2003-fixb",info="static 2003 mixture model on group from akm")
  rkiv0.put(rs,list(res_mixt_akm=res_mixt_akm,res_mixt_blm=res_mixt_blm,grps_kmean=grps_kmean,grp_akm=grp_akm))

  # trying the GIBBS sampling
  sim    = grouping.append(sim,grps_kmean$best_cluster)
  pl_blm = m2.proba.getlcasspr(sim,10)$class_pr
  sim    = grouping.append(sim,grp_akm)
  pl_akm = m2.proba.getlcasspr(sim,10)$class_pr

  rr.gibbs.akm = m2.proba.gibbs(sim,res_mixt_akm$model,pl_akm,nfirms_to_update = 1,grp_start = grp_akm,maxiter = 1000)
  rr.gibbs.blm = m2.proba.gibbs(sim,res_mixt_blm$model,pl_blm,nfirms_to_update = 1,grp_start = grps_kmean$best_cluster,maxiter = 1000)

  rr.akm = data.table(rr.gibbs.akm$rr)
  rr.akm$model="akm"
  rr.blm = data.table(rr.gibbs.blm$rr)
  rr.blm$model="blm"
  rr.gibbs = rbind(rr.akm,rr.blm)

  rr.gibbs[,lik:= likm + liks][,step:=NULL]
  rrm = data.table(melt(rr.gibbs,c("i","model")) )

  ggplot(rrm,aes(x=i,y=value,color=factor(model))) +
    geom_line() + facet_wrap(~variable, scales = "free") + theme_bw()


  # use last classification of AKM, compare to BLM classification
  # or use the BLM clusters
  cl = makeCluster(17)
  clusterEvalQ(cl,require(blm))

  #sim    = grouping.append(sim,rr.gibbs.akm$last_grp)
  #res_mixt_blm_at_akm_last_gibbs = m2.mixt.estimate.all(sim,nk=6,ctrl,cl)
  sim    = grouping.append(sim,grps_kmean$best_cluster)
  sim$sdata[,y1:=y1_bu]
  iter_res = res_mixt_blm
  iter_pl  = m2.proba.getlcasspr(sim,10)$class_pr
  iter_grp = grps_kmean$best_cluster

  # look re-calssifcation/estimation
  rr.all  = data.frame()
  rr.all2 = data.frame()
  for (i in 177:200) {
    # run gibbs 200 firms
    iter_gibbs = m2.proba.gibbs(sim,iter_res$model,iter_pl,nfirms_to_update = 1,grp_start = iter_grp,maxiter = 200)
    iter_grp   = iter_gibbs$last_grp
    # apply latest classification
    sim        = grouping.append(sim,iter_gibbs$last_grp)
    iter_pl    = m2.proba.getlcasspr(sim,10)$class_pr
    # estimate BLM

    cl = makeCluster(10)
    clusterEvalQ(cl,require(blm))
    iter_res = m2.mixt.estimate.all(sim,nk=6,ctrl,cl)
    stopCluster(cl)

    rr.all = rbind(rr.all,data.frame(iter_res$vdec$stats))
    iter_gibbs$rr$outloop=i
    rr.all2 = rbind(rr.all2,iter_gibbs$rr)
  }


  grp2 = rr.gibbs.akm$last_grp
  grp1 = grps_kmean$best_cluster

  save.image(file="../../../data/workspace_proba_approach_akm.Rdata")
  save.image(file="../../../data/workspace_proba_approach_start_blm.Rdata")

  rr.gibbs.blm2 = m2.proba.gibbs(sim,res_mixt_blm$model,pl_blm,nfirms_to_update = 1,grp_start = grps_kmean$best_cluster,maxiter = 1000,mobility_smooth = 0.1)


  rr.gibbs.akm = m2.proba.gibbs(sim,res_mixt_akm$model,pl_akm,nfirms_to_update = 1,grp_start = grp_akm,maxiter = 1000,mobility_smooth = )

  load("../../../data/workspace_proba_approach_akm.Rdata")
  gibbs.all = list()
  gibbs.all$akm$vdec = rr.all
  gibbs.all$akm$liks = rr.all2
  load("../../../data/workspace_proba_approach_start_blm.Rdata")
  gibbs.all$blm$vdec = rr.all
  gibbs.all$blm$liks = rr.all2

  rs    = rkiv0.start("m2-proba-gibbs-d2003-res",info="gibbs estimation results starting from BLM and AKM")
  rkiv0.put(rs,gibbs.all)

  res.all$blm.lik$model="blm"
  res.all$akm.lik$model="akm"
  res.all$blm.lik[,t:=1:.N]
  res.all$akm.lik[,t:=1:.N]

  rr = rbind(res.all$blm.lik,res.all$akm.lik)
  rrm = data.table(melt(rr,c("model","i","step","updates","t")))
  ggplot(rrm,aes(x=t,y=value,color=model)) +
    facet_wrap(~variable,scale="free") +geom_line() + theme_bw()

}


# ========== FIXED EFFECT MODEL ===========

#' This checks the effect of discretizing heterogeneity
#' on the AKM decomposition
server.fe.effect_of_discrete <- function() {
  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]
  sdata=sdata[move==FALSE]
  sim = list(jdata=jdata,sdata=sdata)
  rm(sdata,jdata)

  # ru AKM on jdata to extract firm FE
  # firm.fe = m2.fe.firms(sim$jdata)
  res0 = m2.trace.estimate(sim)

  # classify according to AKM
  require(Ckmeans.1d.dp)
  adata = rbind(sim$sdata[,list(wid,y1,f1,y2,f2)],sim$jdata[,list(wid,y1,f1,y2,f2)])

  setkey(res0$fids,f1)
  setkey(adata,f1)
  adata[, psi1 := res0$fids[adata,psi]]
  setkey(adata,f2)
  adata[, psi2 := res0$fids[adata,psi]]
  adata[, alpha := 0.5*(y1-psi1 + y2-psi2)]

  adata = adata[!is.na(psi1*alpha)]
  dd = adata[, list(nf=Inf,nk=Inf,vpsi = var(psi1) , valpha= var(alpha), cov = cov(alpha,psi1))]

  rr = data.frame()
  rr = rbind(rr,dd)

  for (nf in c(3,10,50)) {

    setkey(adata,wid)

    # discretize psi1
    clusters       = Ckmeans.1d.dp(adata$psi1, nf)
    adata$psi_d    = clusters$centers[clusters$cluster]

    for (nk in c(3,6,10)) {
      # discretize alpha
      clusters       = Ckmeans.1d.dp(adata$alpha, nk)
      adata$alpha_d  = clusters$centers[clusters$cluster]

      dd = adata[, data.frame(nf=nf,nk=nk,vpsi = var(psi_d) , valpha= var(alpha_d), cov = cov(alpha_d,psi_d))]
      rr = rbind(rr,dd)
    }
  }

  rr = data.table(rr)
  rs = rkiv0.start("m2-akm-effect-of-discretization",info="checking how different terms change when discretizing")
  rkiv0.put(rs,rr)

  setkey(rr,nk,nf)
  # secondly, we can estimate
}


sever.fe.trace <- function() {

  load(sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))
  sdata[,x:=1][,y1_bu:=y1]
  sdata=sdata[move==FALSE]
  sim = list(jdata=jdata,sdata=sdata)
  rm(sdata,jdata)

  nmseq = c(100,150,Inf)
  res = m2.trace.blmhyb(sim,use_con = TRUE)

  rs = rkiv0.start("m2-blmhybrid",info="running BLM/AKM while giving large firms their won cluster")
  rkiv0.put(rs,res$rr)

}





