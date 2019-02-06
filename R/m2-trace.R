#' some work on trace formula

#' Extracts the largest connected set from data using f1,f2
#' movers
#' @export
get.largest.conset.fid <- function(jdata) {
  # combine th 2 directions
  jdata2 = rBind(jdata[,list(f1,f2)],jdata[,list(f1=f2,f2=f1)])
  AD = pmin(acast(jdata2[,.N,list(f1,f2)],value.var = "N",f1~f2,drop=FALSE,fill=0),1)
  # compute connected sets
  cs = conComp(AD,2)
  # extract largest
  cs = data.table(f1=names(cs),set=as.numeric(cs))
  cs.size = cs[,.N,set][order(-N)]
  return(cs[  set == cs.size$set[1]  , f1])
}

#' Extracts the largest connected set from data using f1,f2
#' movers
#' @export
get.largest.conset.clus <- function(jdata) {
  # combine th 2 directions
  jdata2 = rBind(jdata[,list(j1,j2)],jdata[,list(j1=j2,j2=j1)])
  AD = pmin(acast(jdata2[,.N,list(j1,j2)],value.var = "N",j1~j2,drop=FALSE,fill=0),1)
  # compute connected sets
  cs = conComp(AD,2)
  # extract largest
  cs = data.table(j1=as.integer(names(cs)),set=as.numeric(cs))
  cs.size = cs[,.N,set][order(-N)]
  return(cs[  set == cs.size$set[1]  , j1])
}


#' Extracts the largest leave-out connected set from data using f1,f2
#' movers
#' @export
get.largest.leaveoutset.fid <- function(jdata) {

  # we loop multiple times
  for (rep in 1:20) {

    # get the connected set
    f1s   = get.largest.conset.fid(jdata)
    jdata = jdata[f1%in%f1s][f2%in%f1s]

    # remove firms with 1 mover
    for (i in 1:10) {
      f0s   = jdata[,list(f1=c(f1,f2),.N)][,.N,f1][N==1,f1]
      if (length(f0s)==0)  break;
      jdata = jdata[!f1%in%f0s][!f2%in%f0s]
    }

    # extract articulation firms
    G   = graph(c(jdata[,rbind(f1,f2)]),directed = F)
    L   = articulation.points(G)
    L   = names(V(G))[L]

    # for each articulation firm, check removing movers
    bad_movers = c()
    for (fi in L) {
      # find edges for this firm
      II1 = jdata[,list(1:.N,f1)][f1==fi,V1]
      II2 = jdata[,list(1:.N,f2)][f2==fi,V1]
      II = union(II1,II2)
      II = setdiff(II,bad_movers)
      for (i in II) {
        Gsub = delete.edges(G, i)
        ng = length( decompose.graph(Gsub) )
        if (ng>1) bad_movers = c(bad_movers,i);
      }
    }

    if (length(bad_movers)==0) break;
    jdata = jdata[setdiff(1:.N,bad_movers)]
  }

  return(jdata)
}

#' Extracts the largest leave-out connected set from data using f1,f2
#' movers
#' @export
get.largest.leaveoutset.clus <- function(jdata) {

  # we loop multiple times
  for (rep in 1:20) {

    # get the connected set
    j1s   = get.largest.conset.clus(jdata)
    jdata = jdata[j1%in%j1s][j2%in%j1s]

    # remove firms with 1 mover
    for (i in 1:10) {
      j0s   = jdata[,list(j1=c(j1,j2),.N)][,.N,j1][N==1,j1]
      if (length(j0s)==0)  break;
      jdata = jdata[!j1%in%j0s][!j2%in%j0s]
    }

    # extract articulation firms
    G   = graph(c(jdata[,rbind(j1,j2)]),directed = F)
    L   = articulation.points(G)
    L   = names(V(G))[L]

    # for each articulation firm, check removing movers
    bad_movers = c()
    for (ji in L) {
      # find edges for this firm
      II1 = jdata[,list(1:.N,j1)][j1==ji,V1]
      II2 = jdata[,list(1:.N,j2)][j2==ji,V1]
      II = union(II1,II2)
      II = setdiff(II,bad_movers)
      for (i in II) {
        Gsub = delete.edges(G, i)
        ng = length( decompose.graph(Gsub) )
        if (ng>1) bad_movers = c(bad_movers,i);
      }
    }

    if (length(bad_movers)==0) break;
    jdata = jdata[setdiff(1:.N,bad_movers)]
  }

  return(jdata)
}



#' create a model for testing trace estimation
#' @export
m2.trace.new <- function(nf = 200 , nm = 10000, eps_sd = 1.5) {
  p = list(dsize = c(exp(3.199292), 1+exp(-1.247662)) , nf = nf , nm = nm, eps_sd = eps_sd  )
  S = rpareto(p$nf,p$dsize[1],p$dsize[2])
  psi = rnorm(p$nf)

  model     = copy(p)
  model$S   = S
  model$psi = psi

  return(model)
}

#' simulates for trace estimation
#' @export
m2.trace.simulate.old <- function(model) {
  JJ = array(0,c(model$nm,model$nf))
  AA = array(0,c(model$nf,model$nf))
  F1 = rep(0,model$nm)
  F2 = rep(0,model$nm)
  for (i in 1:model$nm) {
    ii = sample.int(model$nf,2,prob = model$S)
    JJ[i,ii[1]] = 1
    JJ[i,ii[2]] = -1
    AA[ii[1],ii[2]]=1
    AA[ii[2],ii[1]]=1
    F1[i] = ii[1]
    F2[i] = ii[2]
  }
  colnames(AA) = 1:model$nf
  rownames(AA) = 1:model$nf

  D = JJ %*% model$psi + rnorm(model$nm)*model$eps_sd

  return(list(JJ=JJ,AA=AA,D=D,F1=F1,F2=F2))
}

#' simulates for trace estimation
#' @export
m2.trace.simulate <- function(model) {
  F1   = rep(0,model$nm)
  F2   = rep(0,model$nm)
  psi1 = rep(0,model$nm)
  psi2 = rep(0,model$nm)
  for (i in 1:model$nm) {
    ii = sample.int(model$nf,2,prob = model$S)
    F1[i] = ii[1]
    F2[i] = ii[2]
    psi1[i] = model$psi[ii[1]]
    psi2[i] = model$psi[ii[2]]
  }

  jdata = data.table(f1=paste(F1),f2=paste(F2),psi1=psi1,psi2=psi2)
  jdata[, y1 := psi1 + rnorm(.N)*model$eps_sd/sqrt(2)]
  jdata[, y2 := psi2 + rnorm(.N)*model$eps_sd/sqrt(2)]

  sdata = data.table(f1=paste(1:model$nf),psi1=model$psi,size=model$S)
  sdata = sdata[,list(y1 = psi1 + rnorm(ceiling(size))*model$eps_sd/sqrt(2),
                      y2 = psi1 + rnorm(ceiling(size))*model$eps_sd/sqrt(2),f2=f1), list(f1,size,psi1) ]
  sdata$x=1

  return(list(sdata=sdata,jdata=jdata))
}



#' gets the connected set, then
#' @export
m2.trace.estimate <- function(sim, model0=NA,hetero=FALSE,within_re=FALSE) {

  stats = list()
  stats$hetero = hetero
  stats$total_number_of_firms   = sim$jdata[,length(unique(f1,f2))]
  stats$total_number_of_stayers = sim$sdata[,.N]
  stats$total_number_of_movers  = sim$jdata[,.N]
  stats$total_logwage_var       = var(c(sim$sdata$y1,sim$jdata$y1))
  stats$total_btw_firm          = sim$sdata[,list(mean(y1),.N),f1][,wt.var(V1,N)]

  # EXTRACT CONNECTED SET
  jdata = sim$jdata
  if (hetero==F) {
    f1s   = get.largest.conset.fid(jdata)
    flog.info("connected set %i/%i",length(f1s),length(unique(sim$sdata$f1)))
    jdata = jdata[f1 %in% f1s][f2 %in% f1s]
  } else {
    jdata = get.largest.leaveoutset.fid(jdata)
    f1s = jdata[,unique(c(f1,f2))]
    flog.info("leave-out connected set %i/%i",length(f1s),length(unique(sim$sdata$f1)))
  }

  # index firms with integers
  fids = data.table(f1=f1s,nfid=1:length(f1s))
  setkey(fids,f1)
  setkey(jdata,f1)
  jdata[,f1i := fids[jdata,nfid]]
  setkey(jdata,f2)
  jdata[,f2i := fids[jdata,nfid]]

  # extract size in stayers, and mean
  fsize = rbind(sim$sdata[,list(y1,f1)],sim$jdata[,list(y1,f1)])[,list(N=.N,mw=mean(y1)),f1]
  setkey(fids,f1)
  setkey(fsize,f1)
  fids[,size := fsize[fids,N]]
  fids[,mw   := fsize[fids,mw]]
  fids[is.na(size),size:=0] # pads missing size with 0
  fids[size==0, mw:=0] # pads missing size with 0

  # CONSTRUCT SPARSE DESIGN MATRIX + NORMALIZATION
  nf = pmax(max(jdata$f1i),max(jdata$f2i))
  dd = jdata[,list(m1 = y2 - y1,f1i,f2i)]
  dd = rbind(dd,data.frame(m1=0,f1i=1,f2i=1))
  N = nrow(dd)
  dd[,v1:=-1][,v2:=1]
  dd[,c1:=f1i][,c2:=f2i]
  dd[N,v1:=0]
  JJ = sparseMatrix2(1:N,dd$c1,dd$v1,c(N,nf)) + sparseMatrix2(1:N,dd$c2,dd$v2,c(N,nf))
  S  = fids[order(nfid),size]

  stats$set_number_of_firms   = nf
  stats$set_number_of_stayers = sum(S)-(N-1)
  stats$set_number_of_movers  = N-1
  stats$set_logwage_var       = var(c(sim$sdata[f1 %in% fids$f1,y1],sim$jdata$y1))
  stats$total_btw_firm        = sim$sdata[f1 %in% fids$f1,list(mean(y1),.N),f1][,wt.var(V1,N)]

  # COMPUTE INVERSE OF DESIGN MATRIX
  M = SparseM::t(JJ) %*% JJ
  Minv = SparseM::solve(M,nnzlmax=1e8,tmpmax=1e8,nsubmax=1e8)

  # compute firms FE
  psi = as.numeric(SparseM::as.matrix(Minv %*% ( SparseM::t(JJ) %*% dd$m1)))
  E   = dd$m1 - JJ%*%psi
  E   = E[1:(N-1)]

  # extract random effect variance
  if (within_re==TRUE) {
    omega = m2.trace.reExt(jdata,psi)
  }

  if (!any(is.na(model0))) {
    psi0 = model0$psi[as.integer(fids[order(nfid),f1])]
    psi0 = psi0-psi0[1]
    flog.info("corr= %f", cor(psi0,psi))
  }

  # extract homoskedastic error
  var_e = var(E)*nrow(sim$jdata)/(nrow(sim$jdata)-nf)
  stats$error_var_homo = var_e

  # COMPUTE THE TRACE FORMULA - USING THE WEIGHTING OF STAYERS
  tr_correction = var_e*( sum( SparseM::diag(Minv)*S )/sum(S) - sum( S* (Minv %*% rep(1/nf,nf)) )/(sum(S))  )
  stats$trace_term_homo = tr_correction
  stats$trace_term_hetero = NA
  stats$error_var_hetero = NA

  # heteroskedastic variance using BLM groups
  if (hetero==TRUE) {

    S2 = S/sum(S)
    V  = Minv %*% SparseM::t(JJ)
    # we construct a slightly different trace using esitmate of individual variance
    # suing KKS.
    P       = sColSums(SparseM::t(JJ) * V)
    S_i     = as.numeric(dd$m1 * (dd$m1 - JJ%*%psi)/(1-P))
    stats$error_var_hetero = mean(S_i)

    # Finally we need to compute the full correction
    B = as.numeric(sColSums(V *  diag(S2) %*% V) - sColSums(diag(S2) %*% V ) ^2)
    tr_correction_hetero = sum( S_i * B)

    flog.info("tr0=%f tr1=%f var0=%f var1=%f",tr_correction, tr_correction_hetero,var_e,mean(S_i))
    tr_correction = tr_correction_hetero
    stats$trace_term_hetero = tr_correction_hetero
  }

  stats$trace_term = tr_correction

  # COMPUTE THE VARIANCE OF PSI
  fids[, psi := psi[nfid]]
  if (within_re==TRUE) {
    fids[,omega := omega[nfid]]
  }
  var_psi_hat = fids[,wt.var(psi,size)]
  tot = sim$sdata[,var(y1)]
  btw = sim$sdata[,list(mean(y1),.N),f1][,wt.var(V1,N)]

  stats$psi_var = fids[,wt.var(psi,size)]

  if (!any(is.na(model0))) {
    fids[, psi0 := psi0[nfid]]
    rm(psi0)
    flog.info("var_true=%f  var_akm=%f var2=%f trace=%f ", fids[,wt.var(psi0,size)], fids[,wt.var(psi,size)], fids[,wt.var(psi,size)]- tr_correction,tr_correction)
  } else {
    flog.info("tot=%f btwf=%f var_akm=%f var2=%f trace=%f ",tot,btw, fids[,wt.var(psi,size)], fids[,wt.var(psi,size)]- tr_correction,tr_correction)
  }

  res = list(fids=fids,eps_sd = sqrt(var_e), var_psi= var_psi_hat, stats=stats)
}

#' Estimates a two-sided heterogeneity model at the
#' cluster level using (j1,j2).
#' @param hetero if TRUE, extracts the leave-on-out connected set
#' @param within_re whether to compute the within cluster random firm heterogeneity
#' @export
m2.trace.estimate.clus <- function(sim, model0=NA,hetero=FALSE,within_re=FALSE) {

  stats = list()
  stats$hetero = hetero
  stats$total_number_of_firms      = sim$jdata[,length(unique(f1,f2))]
  stats$total_number_of_clusters   = sim$jdata[,length(unique(j1,j2))]
  stats$total_number_of_stayers    = sim$sdata[,.N]
  stats$total_number_of_movers     = sim$jdata[,.N]
  stats$total_logwage_var          = var(c(sim$sdata$y1,sim$jdata$y1))
  stats$total_btw_firm             = sim$sdata[,list(mean(y1),.N),f1][,wt.var(V1,N)]

  # EXTRACT CONNECTED SET
  jdata = sim$jdata
  if (hetero==F) {
    j1s   = get.largest.conset.clus(jdata)
    flog.info("connected set %i/%i clusters",length(j1s),length(unique(sim$sdata$j1)))
    jdata = jdata[j1 %in% j1s][j2 %in% j1s]
  } else {
    jdata = get.largest.leaveoutset.clus(jdata)
    j1s = jdata[,unique(c(j1,j2))]
    flog.info("leave-out connected set %i/%i",length(j1s),length(unique(sim$sdata$j1)))
  }

  # index groups with integers (should be the case already)
  jids = data.table(j1=j1s,nfid=1:length(j1s))
  setkey(jids,j1)
  setkey(jdata,j1)
  jdata[,j1i := jids[jdata,nfid]]
  setkey(jdata,j2)
  jdata[,j2i := jids[jdata,nfid]]

  # extract size in stayers, and mean
  fsize = rbind(sim$sdata[,list(y1,j1)],sim$jdata[,list(y1,j1)])[,list(N=.N,mw=mean(y1)),j1]
  setkey(jids,j1)
  setkey(fsize,j1)
  jids[,size := fsize[jids,N]]
  jids[,mw   := fsize[jids,mw]]
  jids[is.na(size),size:=0] # pads missing size with 0
  jids[size==0, mw:=0] # pads missing size with 0

  # CONSTRUCT SPARSE DESIGN MATRIX + NORMALIZATION
  nf = pmax(max(jdata$j1i),max(jdata$j2i))
  dd = jdata[,list(m1 = y2 - y1,j1i,j2i)]
  dd = rbind(dd,data.frame(m1=0,j1i=1,j2i=1))
  N = nrow(dd)
  dd[,v1:=-1][,v2:=1]
  dd[,c1:=j1i][,c2:=j2i]
  dd[N,v1:=0]
  JJ = sparseMatrix2(1:N,dd$c1,dd$v1,c(N,nf)) + sparseMatrix2(1:N,dd$c2,dd$v2,c(N,nf))
  S  = jids[order(nfid),size]

  stats$set_number_of_firms   = nf
  stats$set_number_of_stayers = sum(S)-(N-1)
  stats$set_number_of_movers  = N-1
  stats$set_logwage_var       = var(c(sim$sdata[j1 %in% jids$j1,y1],sim$jdata$y1))
  stats$total_btw_firm        = sim$sdata[j1 %in% jids$j1,list(mean(y1),.N),j1][,wt.var(V1,N)]

  # COMPUTE INVERSE OF DESIGN MATRIX
  M = SparseM::t(JJ) %*% JJ
  Minv = SparseM::solve(M,nnzlmax=1e8,tmpmax=1e8,nsubmax=1e8)

  # compute firms FE
  psi = as.numeric(SparseM::as.matrix(Minv %*% ( SparseM::t(JJ) %*% dd$m1)))
  E   = dd$m1 - JJ%*%psi
  E   = E[1:(N-1)]

  # extract random effect variance
  if (within_re==TRUE) {
    jdata[,psi1 := psi[j1i]]
    jdata[,psi2 := psi[j2i]]
    omega = m2.trace.reExt(jdata)
  }

  if (!any(is.na(model0))) {
    psi0 = model0$psi[as.integer(jids[order(nfid),j1])]
    psi0 = psi0-psi0[1]
    flog.info("corr= %f", cor(psi0,psi))
  }

  # extract homoskedastic error
  var_e = var(E)*nrow(sim$jdata)/(nrow(sim$jdata)-nf)
  stats$error_var_homo = var_e

  # COMPUTE THE TRACE FORMULA - USING THE WEIGHTING OF STAYERS
  tr_correction = var_e*( sum( SparseM::diag(Minv)*S )/sum(S) - sum( S* (Minv %*% rep(1/nf,nf)) )/(sum(S))  )
  stats$trace_term_homo = tr_correction
  stats$trace_term_hetero = NA
  stats$error_var_hetero = NA

  # heteroskedastic variance using BLM groups
  if (hetero==TRUE) {

    S2 = S/sum(S)
    V  = Minv %*% SparseM::t(JJ)
    # we construct a slightly different trace using esitmate of individual variance
    # suing KKS.
    P       = sColSums(SparseM::t(JJ) * V)
    S_i     = as.numeric(dd$m1 * (dd$m1 - JJ%*%psi)/(1-P))
    stats$error_var_hetero = mean(S_i)

    # Finally we need to compute the full correction
    B = as.numeric(sColSums(V *  diag(S2) %*% V) - sColSums(diag(S2) %*% V ) ^2)
    tr_correction_hetero = sum( S_i * B)

    flog.info("tr0=%f tr1=%f var0=%f var1=%f",tr_correction, tr_correction_hetero,var_e,mean(S_i))
    tr_correction = tr_correction_hetero
    stats$trace_term_hetero = tr_correction_hetero
  }

  stats$trace_term = tr_correction

  # COMPUTE THE VARIANCE OF PSI
  jids[, psi := psi[nfid]]
  if (within_re==TRUE) {
    jids[,omega := omega[nfid]]
  }
  var_psi_hat = jids[,wt.var(psi,size)]
  tot = sim$sdata[,var(y1)]
  btw = sim$sdata[,list(mean(y1),.N),j1][,wt.var(V1,N)]

  stats$psi_var = jids[,wt.var(psi,size)]

  if (!any(is.na(model0))) {
    jids[, psi0 := psi0[nfid]]
    rm(psi0)
    flog.info("var_true=%f  var_akm=%f var2=%f trace=%f ", jids[,wt.var(psi0,size)], jids[,wt.var(psi,size)], jids[,wt.var(psi,size)]- tr_correction,tr_correction)
  } else {
    flog.info("tot=%f btwf=%f var_akm=%f var2=%f trace=%f ",tot,btw, jids[,wt.var(psi,size)], jids[,wt.var(psi,size)]- tr_correction,tr_correction)
  }

  res = list(jids=jids,eps_sd = sqrt(var_e), var_psi= var_psi_hat, stats=stats)
}


check.data <- function(sim) {
  sf1 = sim$sdata[,unique(f1)]
  jf1 = sim$jdata[,unique(f1)]
  jf2 = sim$jdata[,unique(f2)]

  # compute firms in sdata but not in jdata
  flog.info(" %i firms in sdata but not in jdata",   length(setdiff( sf1, c(jf1,jf2) )  ))
  flog.info(" %i firms in jdata.f1 but not in sdata",length(setdiff( jf1, sf1 )))
  flog.info(" %i firms in jdata.f2 but not in sdata",length(setdiff( jf2, sf1 )))
  flog.info(" %i firms in p2 but not in p1",length(setdiff( jf2, c(sf1,jf1 ))))
}


#' estimates a hybrid between BLM and AKM
#'
#' the function returns a bunch of stuff. It returns info about different sets:
#'  - set0 is just the data it starts with at each nm
#'  - set1 is AKM connected set
#'  - set2 is the leave-out connected set
#'  - set3 is the connected set only among firms with at least nm movers
#'
#' @export
m2.trace.blmhyb <- function(sim,use_con=FALSE,nclus=10,nm_list=c(seq(1,10),seq(20,50,by=5),100,150,Inf)) {

  # extract connected set
  if (use_con==T) {
    f0s   = get.largest.conset.fid(sim$jdata)
    jdata = sim$jdata[f1%in%f0s][f2%in%f0s]
    sdata = sim$sdata[f1%in%f0s]
  }

  # save a copy of the data
  sim.copy = copy(sim)

  # next we group differently,
  rr2 = data.frame()
  for (nm in nm_list) {

    stats = list()
    sim = copy(sim.copy)
    stats$nm = nm

    # select firms from AKM
    if (use_con==T) {
      sim$sdata = sim$sdata[f1 %in% f0s]
      sim$jdata = sim$jdata[f1 %in% f0s][f2 %in%f0s]
    }

    # ------- (0) STATS ON CONNECTED SET ------ #
    stats$set0_nf       = length(unique(c(sim$jdata[,unique(c(f1,f2))],sim$sdata$f1)))
    stats$set0_nw       = sim$sdata[,.N] + sim$jdata[,.N]
    stats$set0_wage_var = sim$sdata[,var(y1)]
    stats$set0_btw_firm = sim$sdata[,list(mean(y1),.N),f1][,wt.var(V1,N)]

    # we find firms with at least nm movers
    large.firms               = sim$jdata[,list(f1=c(f1,f2))][,.N,f1][N>=nm,f1]
    stats$firms_with_own_type = length(large.firms)

    # we cluster the other firms
    if (nm>1) {
      ms    = grouping.getMeasures(list(sdata = sim$sdata[!(f1 %in% large.firms)],jdata=sim$jdata),"ecdf",Nw=20,y_var = "y1")
      grps  = grouping.classify.once(ms,k = nclus,nstart = 1000,iter.max = 200,step=100)

      # we create the full classification
      clus0 = grps$best_cluster
      cluster_with_many_firms = sort(unique(clus0))
      clus1 = 1:length(large.firms) + nclus
      names(clus1) = large.firms
      clus = c(clus0,clus1)
    } else {
      cluster_with_many_firms = c()
      clus = 1:length(large.firms)
      names(clus) = large.firms
    }

    # append the clusters
    sim  = grouping.append(sim,clus,drop=T,sort=FALSE)

    # ------- (1) AKM + ANDREWS ------ #
    fids_con   = get.largest.conset.clus(sim$jdata)
    sim$jdata  = sim$jdata[j1%in%fids_con][j2%in%fids_con]
    sim$sdata  = sim$sdata[j1%in%fids_con][j2%in%fids_con]

    stats$set1_nf = length(unique(c(sim$jdata[,unique(c(f1,f2))],sim$sdata$f1)))
    stats$set1_nw = sim$sdata[,.N] + sim$jdata[,.N]
    stats$set1_wage_var = sim$sdata[,var(y1)]
    stats$set1_btw_firm = sim$sdata[,list(mean(y1),.N),f1][,wt.var(V1,N)]

    res_hybrid    = m2.trace.estimate.clus(sim,within_re = TRUE)

    stats$set1_psi_var   = res_hybrid$var_psi
    stats$set1_cov       = res_hybrid$jids[,wt.cov(mw-psi,psi,size)]
    stats$set1_trace     = res_hybrid$stats$trace_term
    stats$set1_omega_var = res_hybrid$jids[j1 %in% cluster_with_many_firms,wt.mean(omega,size)]
    stats$set1_omega_var_pos = res_hybrid$jids[j1 %in% cluster_with_many_firms,wt.mean(pmax(omega,0),size)]

    if (nm>1) {
      stats$set1_psi_var_own  = res_hybrid$jids[!(j1 %in% cluster_with_many_firms),wt.var(psi,size)]
    } else {
      stats$set1_psi_var_own = stats$set1_psi_var
    }

    # ------- (2) HETERO TRACE ------ #
    sim$jdata = get.largest.leaveoutset.clus(sim$jdata)
    jids_con  = sim$jdata[,unique(c(j1,j2))]
    sim$sdata = sim$sdata[j1%in%jids_con][j2%in%jids_con]

    stats$set2_nf = length(unique(c(sim$jdata[,unique(c(f1,f2))],sim$sdata$f1)))
    stats$set2_nw = sim$sdata[,.N] + sim$jdata[,.N]
    stats$set2_wage_var = sim$sdata[,var(y1)]
    stats$set2_btw_firm = sim$sdata[,list(mean(y1),.N),f1][,wt.var(V1,N)]

    res_hybrid    = m2.trace.estimate.clus(sim,hetero = T,within_re = TRUE)
    res_hybrid$stats$nm = nm
    stats$set2_psi_var  = res_hybrid$var_psi
    stats$set2_cov      = res_hybrid$jids[,wt.cov(mw-psi,psi,size)]
    stats$set2_trace    = res_hybrid$stats$trace_term

    # ------- (3) AKM ON FIRM WITH OWN TYPE ------ #
    sim = copy(sim.copy)
    sim$sdata = sim$sdata[f1 %in% large.firms]
    sim$jdata = sim$jdata[f1 %in% large.firms][f2 %in% large.firms]

    stats$set3_nf       = length(unique(c(sim$jdata[,unique(c(f1,f2))],sim$sdata$f1)))
    stats$set3_nw       = sim$sdata[,.N] + sim$jdata[,.N]
    stats$set3_wage_var = sim$sdata[,var(y1)]
    stats$set3_btw_firm = sim$sdata[,list(mean(y1),.N),f1][,wt.var(V1,N)]

    stats$set3_psi_var = NA
    stats$set3_cov = NA
    stats$set3_trace = NA

    if (sim$jdata[,length(unique(c(f1,f2)))]>1) {
      fids_con   = get.largest.conset.fid(sim$jdata)
      sim$jdata  = sim$jdata[f1%in%fids_con][f2%in%fids_con]
      sim$sdata  = sim$sdata[f1%in%fids_con][f2%in%fids_con]

      stats$set3_nf = length(unique(c(sim$jdata[,unique(c(f1,f2))],sim$sdata$f1)))
      stats$set3_nw = sim$sdata[,.N] + sim$jdata[,.N]
      stats$set3_wage_var = sim$sdata[,var(y1)]

      res_hybrid  = m2.trace.estimate(sim)
      res_hybrid$stats$nm=nm
      stats$set3_psi_var = res_hybrid$var_psi
      stats$set3_cov     = res_hybrid$fids[,wt.cov(mw-psi,psi,size)]
      stats$set3_trace   = res_hybrid$stats$trace_term
    }

    # ------- (4) AKM-HETERO ON FIRM WITH OWN TYPE ------ #
    sim = copy(sim.copy)
    sim$sdata = sim$sdata[f1 %in% large.firms]
    sim$jdata = sim$jdata[f1 %in% large.firms][f2 %in% large.firms]

    stats$set4_nf       = length(unique(c(sim$jdata[,unique(c(f1,f2))],sim$sdata$f1)))
    stats$set4_nw       = sim$sdata[,.N] + sim$jdata[,.N]
    stats$set4_wage_var = sim$sdata[,var(y1)]
    stats$set4_btw_firm = sim$sdata[,list(mean(y1),.N),f1][,wt.var(V1,N)]

    stats$set4_psi_var = NA
    stats$set4_cov = NA
    stats$set4_trace = NA

    if (sim$jdata[,length(unique(c(f1,f2)))]>1) {
      # leave-one-out
      sim$jdata = get.largest.leaveoutset.clus(sim$jdata)
      jids_con  = sim$jdata[,unique(c(j1,j2))]
      sim$sdata = sim$sdata[j1%in%jids_con][j2%in%jids_con]

      stats$set4_nf = length(unique(c(sim$jdata[,unique(c(f1,f2))],sim$sdata$f1)))
      stats$set4_nw = sim$sdata[,.N] + sim$jdata[,.N]
      stats$set4_wage_var = sim$sdata[,var(y1)]

      res_hybrid  = m2.trace.estimate(sim,hetero = T)
      res_hybrid$stats$nm=nm
      stats$set4_psi_var = res_hybrid$var_psi
      stats$set4_cov     = res_hybrid$fids[,wt.cov(mw-psi,psi,size)]
      stats$set4_trace   = res_hybrid$stats$trace_term
    }


    rr2 = rbind(rr2,data.frame(stats))
  }

  return(list(rr=rr2,rr_all=rr))
}

#' Random effect within cluster variance
#' @export
m2.trace.reExt <- function(jdata,subsample=0) {

  nf = max(c(jdata$j1,jdata$j2))

  omega = rep(0,nf)
  for (l1 in 1:nf) {
    # extract the firms with 2 movers in this group
    num= 0; den=0;
    firm_in_grp = jdata[l1==j1,.N,f1][N>=2,f1]
    #flog.info("found %i firms for j1=%i",length(firm_in_grp),l1)
    for (firm_cur in firm_in_grp) {
      dt = jdata[f1==firm_cur]
      if (subsample>0) dt = dt[sample.int(.N,pmin(.N,subsample))];

      # create all pairs, make sure they move to different firms
      F2 = spread(dt$f2,2,nrow(dt))
      YY = spread(dt$y2 - dt$psi2 - dt$y1+dt$psi1,2,nrow(dt))

      # add this to the measure of variance!
      num = num + sum((F2 != t(F2)) * YY * t(YY))/2
      den = den + sum((F2 != t(F2)))/2
    }
    firm_in_grp = jdata[l1==j2,.N,f2][N>=2,f2]
    for (firm_cur in firm_in_grp) {
      dt = jdata[f2==firm_cur]
      if (subsample>0) dt = dt[sample.int(.N,pmin(.N,subsample))];

      # create all pairs, make sure they move to different firms
      F1 = spread(dt$f1,2,nrow(dt))
      YY = spread(dt$y2 - dt$psi2 - dt$y1+dt$psi1,2,nrow(dt))

      # add this to the measure of variance!
      num = num + sum((F1 != t(F1)) * YY * t(YY))/2
      den = den + sum((F1 != t(F1)))/2
    }
    if (den>0) omega[l1] = num/den;
    flog.info("num=%f den=%f var=%f j1=%i",num,den,num/den,l1)
  }

  return(omega)
}





