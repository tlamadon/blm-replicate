#' @export
catf <- function(...) cat(sprintf(...))

#' @export
ssample <- function(dd,s) {
  dd[sample.int(.N,round(.N*s))]
}

#' this is a utility function to generate
#' multidimensional arrays - like the spread function in fortran
#' @export
spread <- function (A, loc, dims) {
  if (!(is.array(A))) {
    A = array(A, dim = c(length(A)))
  }
  adims = dim(A)
  l = length(loc)
  if (max(loc) > length(dim(A)) + l) {
    stop("incorrect dimensions in spread")
  }
  sdim = c(dim(A), dims)
  edim = c()
  oi = 1
  ni = length(dim(A)) + 1
  for (i in c(1:(length(dim(A)) + l))) {
    if (i %in% loc) {
      edim = c(edim, ni)
      ni = ni + 1
    }
    else {
      edim = c(edim, oi)
      oi = oi + 1
    }
  }
  return(aperm(array(A, dim = sdim), edim))
}

hist2 <- function(Y1,Y2,wsup) {
  n = length(wsup)
  H = array(0,c(n,n))
  for (i in 1:n) for (j in 1:n) {
    H[i,j] = sum(  (Y1 < wsup[i]) & (Y2 < wsup[j]) )
  }
  H[n,n] = length(Y1)
  H = H/H[n,n]
  return(H)
}

hist1 <- function(Y1,wsup) {
  n = length(wsup)
  H = array(0,c(n))
  for (i in 1:n) {
    H[i] = sum(  (Y1 < wsup[i]) )
  }
  H[n] = length(Y1)
  H = H/H[n]
  return(H)
}

# allow to use 2 different supports
hist2d <- function(Y1,Y2,wsup) {
  n = length(wsup)
  H = array(0,c(n,n))
  for (i in 1:n) for (j in 1:n) {
    H[i,j] = sum(  (Y1 < wsup[i]) & (Y1 >= wsup[i-1]) & (Y2 < wsup[j]) & (Y1 >= wsup[i-1])  )
  }
  H[n,n] = length(Y1)
  H = H/H[n,n]
  return(H)
}

# smoothed histogram
hist2s <- function(Y1,Y2,wsup,h) {
  n = length(wsup)
  H = array(0,c(n,n))
  for (i in 1:n) for (j in 1:n) {
    H[i,j] = sum(  pnorm( (wsup[i] - Y1 )/h ) *  pnorm( (wsup[j] - Y2 )/h ) )
  }
  H = H/H[n,n]
  return(H)
}

#' @export
rdim <- function(A,...) {
  dd <- list(...);
  if (length(dd)==1) {
    dim(A)<-dd[[1]]
  } else {
    dim(A) <- dd
  }
  return(A)
}



tic.new <- function() {
  t = Sys.time()
  tt = list(all.start=t,last.time=t,loop=0,timers=list())

  tic.toc <- function(name="") {
    t = Sys.time()
    if (name=="") {
      return(tt)
    }

    if (name %in% names(tt$timers)) {
      tm = tt$timers[[name]]
      tm$count = tm$count+1
      tm$total = tm$total + t - tt$last.time
      tt$timers[[name]] = tm
    } else {
      tm = list(count=1,total=t - tt$last.time)
      tt$timers[[name]] = tm
    }
    tt$last.time=t;
    tt <<- tt
  }

  return(tic.toc)
}

#' @export 
addmom <- function(A1,A2,name,rr=data.frame(),type="mean") {
  rr1 = melt(A1)
  rr2 = melt(A2)
  rr = rbind(rr,data.frame(name=name,val1=rr1$value,val2=rr2$val,type=type))
  return(rr)
}

lm.wfitc <- function(XX,YY,rw,C1,C0,meq) {

  S = apply(abs(XX),2,max)
  XX2 = XX*spread(1/S,1,dim(XX)[1])
  C12 = C1*spread(1/S,1,dim(C1)[1])

  XXw      = diag(rw) %*% XX2
  Dq       = t(XXw) %*% XX2
  dq       = t(YY %*% XXw)

  # do quadprod
  fit      = solve.QP(Dq,dq,t(C12),C0,meq)

  # rescale
  fit$solution = fit$solution/S

  return(fit)
}

lm.wfitnn <- function(XX,YY,rw,floor = 0) {

  n = dim(XX)[2]
  XX2 = XX
  #   S = apply(abs(XX),2,max)
  #   XX2 = XX*spread(1/S,1,dim(XX)[1])
  #   C12 = C1*spread(1/S,1,dim(C1)[1])

  XXw      = diag(rw) %*% XX2
  Dq       = t(XXw) %*% XX2
  dq       = t(YY %*% XXw)
  C1       = diag(n)
  C0       = rep(floor,n)

  fit      = qprog(Dq,dq,C1,C0)
  fit$solution = as.numeric(fit$thetahat)

  return(fit)
}

# fits a linear problem with weights under constraints
slm.wfitc <- function(XX,YY,rw,CS,scaling=0) {
  nk = CS$nk
  nf = CS$nf
  YY = as.numeric(YY)
  XX = as.matrix.csr(XX)
  # to make sure the problem is positive semi definite, we add
  # the equality constraints to the XX matrix! nice, no?
  
  if (CS$meq>0) {
    XXb = rbind(XX,  as.matrix.csr(CS$C[1:CS$meq,]))
    YYb = c(YY,CS$H[1:CS$meq])
    rwb  = c(rw,rep(1,CS$meq))
  } else {
    XXb = XX
    YYb = YY
    rwb = rw
  }
  
  t2 = as(dim(XXb)[1],"matrix.diag.csr")
  t2@ra = rwb
  XXw = t2 %*% XXb
  Dq       = SparseM:::as.matrix(SparseM:::t(XXw) %*% XXb)
  dq       = SparseM:::t(YYb %*% XXw)
  
  # scaling
  #
  if (scaling>0) {
    sc <- norm(Dq,"2")^scaling 
  } else {
    sc=1
  }
  
  # do quadprod
  tryCatch({
    fit      = solve.QP(Dq/sc,dq/sc,t(CS$C)/sc,CS$H/sc,CS$meq)
  }, error = function(err) {
    browser()
  })
  
  return(fit)
}



lognormpdf <- function(Y,mu=0,sigma=1)   -0.5 * (  (Y-mu) / sigma )^2   - 0.5 * log(2.0*pi) - log(sigma)  

logsumexp <- function(v) {
  vm = max(v)
  log(sum(exp(v-vm))) + vm
}

mplot <- function(M) {
  mm = melt(M,c('i','j'))
  #mm$i = factor(mm$i)
  #mm$j = factor(mm$j)
  mm = mm[mm$value>0,]
  ggplot(mm,aes(x=j,y=i,fill=value)) + geom_tile() + theme_bw() + scale_y_reverse()
}

mvar <- function(x) { 
  if (length(x)<=1) return(0);
  return(var(x))
}
mcov <- function(x,y) { 
  if (length(x)<=1) return(0);
  return(cov(x,y))
}

#' @export
archive.new <- function(file) {
  desc=data.frame()
  save(desc,file=paste(path_archive,file,sep=""))
}

#' @export
archive.list <- function(file) {
  filename=paste(path_archive,file,sep="")
  # check that
  load(filename)
  return(desc)
}

#' @export
archive.put <- function(...,m,file) {
  params = list(...)
  name=names(params)[[1]]

  filename=paste(path_archive,file,sep="")
  load(filename,verbose=F)

  if (name %in% desc$name) {
    I=which(desc$name==name)
    desc$desc[I] = m
    desc$date[I]= format(Sys.Date(), "%d-%m-%y")
  } else {
    desc = rbind(desc,data.frame(name=paste(name),desc=m,date=format(Sys.Date(), "%d-%m-%y")))
  }
  desc$name=paste(desc$name)
  desc$date = paste(desc$date)
  ll = c('desc',c(desc$name))
  with(params,save(list=ll,file=filename))
}

#' @export
archive.get <- function(name,file) {
  filename=paste(path_archive,file,sep="")
  load(filename)
  
  if (name %in% desc$name) {
    return(get(name))
  } else {
    stop("this variable is not in the archive")
  }
}

#' @export
assert_notna <- function(...) {
  ps = list(...)
  for (n in names(ps)) {
    if (any(is.na(ps[[n]]))) { stop( paste(n,"has NA values"))} 
  }
}

#' @export
wplot <- function(Wm) {
  dd = melt(Wm,c('l','k'))
  ggplot(dd,aes(x=factor(l),color=factor(k),group=factor(k),y=value)) + geom_line() + theme_bw()
}

#' @export
mplot <- function(M) {
  mm = melt(M,c('i','j'))
  #mm$i = factor(mm$i)
  #mm$j = factor(mm$j)
  mm = mm[mm$value>0,]
  ggplot(mm,aes(x=j,y=i,fill=value)) + geom_tile() + theme_bw() + scale_y_reverse()
}

#' @export
pplot <- function(pk0) {
  dd = melt(pk0,c('l','k'))
  ggplot(dd,aes(x=factor(l),y=value,fill=factor(k))) + geom_bar(position="stack",stat = "identity") + theme_bw()
}

# This is adding a moment to compare two models
addmom <- function(A1,A2,name,rr=data.frame(),type="mean") {
  rr1 = melt(A1)
  rr2 = melt(A2)
  rr = rbind(rr,data.frame(name=name,val1=rr1$value,val2=rr2$val,type=type))
  return(rr)
}

#' conditional distribution
#' @export
condquant <- function(Y,X,qs,df=15) {
  XX  = model.matrix(Y ~ bs(X, df=df))
  dd = data.frame()
  for (i in 1:length(qs)) {
    fit = rq(Y ~ bs(X, df=df), tau=qs[i])
    dd = rbind(dd,data.frame(qs=qs[i],y=XX %*% fit$coef,x=X))
  }
  dd
}

condquant2 <- function(Y,X,qs) {
  # using locfit
  
  
  XX  = model.matrix(Y ~ bs(X, df=15))
  dd = data.frame()
  for (i in 1:length(qs)) {
    fit = rq(Y ~ bs(X, df=15), tau=qs[i])
    dd = rbind(dd,data.frame(qs=qs[i],y=XX %*% fit$coef,x=X))
  }
  dd
}

#' reindexed the fids of firms
#' 
#' This does not copy si, it changes it, careful!
#' @export
firm.reindex <- function(sim,fids=c()) {
  # get unique fids
  if (length(fids)==0) {
    fids = unique(c(sim$sdata[,unique(f1)],sim$jdata[,unique(f1)],sim$jdata[,unique(f2)]))  
  } else {
    # we subset
    sim$jdata = sim$jdata[f1 %in% fids][f2 %in% fids]
    sim$sdata = sim$sdata[f1 %in% fids]
  }
  
  fids = data.table(fnew=1:length(fids),f=fids)
  setkey(fids,f)
  setkey(sim$jdata,f1)
  sim$jdata[,f1:=fids[sim$jdata,fnew]]
  setkey(sim$jdata,f2)
  sim$jdata[,f2:=fids[sim$jdata,fnew]]
  setkey(sim$sdata,f1)
  sim$sdata[,f1:=fids[sim$sdata,fnew]]
  
  return(sim)
}

#' order cluster by increasing wage
#' @export
cluster.order <- function(sim) {
  sim$sdata = sim$sdata[!is.na(j1)]
  I = sim$sdata[,mean(y1),j1][order(j1)][,rank(V1)]
  sim$sdata[,j1:=I[j1]][,j2:=I[j2]]
  sim$jdata[,j1:=I[j1]][,j2:=I[j2]]
  return(sim)  

  # sim$sdata = sim$sdata[!is.na(j1)]
  # nf = sim$sdata[,max(j1)]
  # 
  # # dealing with the fact that some cluster might be empty
  # V = rep(-Inf,nf)
  # for (jj in 1:nf) {
  #   V[jj] = sim$sdata[j1==jj,mean(y1)]
  # }
  # I = order(V)
  # 
  # sim$sdata[,j1:=I[j1]][,j2:=I[j2]]
  # sim$jdata[,j1:=I[j1]][,j2:=I[j2]]
  # return(sim)
}
  
#' Compute the LIML estimator in teh case where 
#' the dependent variable is 0.
liml.nodep <- function(X,Z) {
  t = Matrix:::t
  Wz  = t(X) %*% Z %*% solve( t(Z) %*% Z ) %*% t(Z) %*% X
  Wx  = t(X) %*% X
  WW  = solve(Wx) %*% Wz
  dec = eigen(WW)
  b_liml = Re(as.numeric(dec$vectors[,dim(dec$vectors)[2]]))  
  return(b_liml)
}

#' Computes graph connectedness among the movers
#' within each type and returns the smalless value
model.connectiveness <- function(model,all=FALSE) {
  EV = rep(0,model$nk)
  pk1 = rdim(model$pk1,model$nf,model$nf,model$nk)
  dd_post = data.table(melt(pk1,c('j1','j2','k')))
  pp = model$NNm/sum(model$NNm)
  dd_post <- dd_post[, pr_j1j2 := pp[j1,j2],list(j1,j2)  ]
  dd_post <- dd_post[, pr_j1j2k := pr_j1j2*value]
  
  for (kk in 1:model$nk) {
    # compute adjency matrix
    A1 = acast(dd_post[k==kk, list(pr=pr_j1j2k/sum(pr_j1j2k),j2,j1)],j1~j2,value.var = "pr")
    A2 = acast(dd_post[k==kk, list(pr=pr_j1j2k/sum(pr_j1j2k),j2,j1)],j2~j1,value.var = "pr")
    # construct Laplacian
    A = 0.5*A1 + 0.5*A2
    D = diag( rowSums(A)^(-0.5) )
    L = diag(model$nf) - D%*%A%*%D
    EV[kk] = sort(eigen(L)$values)[2]
    #print(eigen(L)$values)
  }
  if (all==TRUE) return(EV);
  return(min(abs(EV)))
}

model.connectiveness2 <- function(model,all=FALSE) {
  EV = rep(0,model$nk)
  pk1 = rdim(model$pk1,10,10,model$nk)
  dd_post = data.table(melt(pk1,c('j1','j2','k')))
  pp = model$NNm/sum(model$NNm)
  dd_post <- dd_post[, pr_j1j2 := pp[j1,j2],list(j1,j2)  ]
  dd_post <- dd_post[, pr_j1j2k := pr_j1j2*value]
  
  for (kk in 1:model$nk) {
    M = acast(dd_post[k==kk, list(pr=pr_j1j2k/sum(pr_j1j2k),j2),j1],j1~j2,value.var = "pr")
    A = acast(dd_post[k==kk, list(pr=pr_j1j2k/sum(pr_j1j2k),j2,j1)],j1~j2,value.var = "pr")
    EV[kk] = -sort(-eigen(M)$values)[2]
  }
  if (all==TRUE) return(EV);
  return(max(abs(EV)))  
}

model.connectiveness3 <- function(model,all=FALSE) {
  EV = rep(0,model$nk)
  pk1 = rdim(model$pk1,10,10,model$nk)
  dd_post = data.table(melt(pk1,c('j1','j2','k')))
  pp = model$NNm/sum(model$NNm)
  dd_post <- dd_post[, pr_j1j2 := pp[j1,j2],list(j1,j2)  ]
  dd_post <- dd_post[, pr_j1j2k := pr_j1j2*value]
  
  for (kk in 1:model$nk) {
    #A = acast(dd_post[k==kk, list(pr=pr_j1j2k/sum(pr_j1j2k),j2,j1)],j1~j2,value.var = "pr")
    M1 = acast(dd_post[k==kk, list(pr=pr_j1j2k/sum(pr_j1j2k),j2),j1],j1~j2,value.var = "pr")
    M2 = acast(dd_post[k==kk, list(pr=pr_j1j2k/sum(pr_j1j2k),j1),j2],j2~j1,value.var = "pr")
    M = 0.5*M1 + 0.5*M2
    EV[kk] = -sort(-eigen(M)$values)[2]
  }
  if (all==TRUE) return(EV);
  return(max(abs(EV)))  
}  

#' joint melt of multiple arrays
multimelt <- function(...) {
  
  # I would want somthing like melt(A,c('n','m')) + melt(A,c('n','m')) + ...
  
  
}

#' compute some stats on data
#' @export
sample.stats <- function(sdata,ydep,type="j1",name=type) {
  sdata[,y_dep := get(ydep)]
  sdata[,jtmp  := get(type)]
  sdata[,y_pred := mean(y_dep),list(k,jtmp)]
  sdata[,y_eps  := y_dep-y_pred,list(k,jtmp)]
  
  # we compute some statistics
  rt = sdata[,list(mean=mean(y_dep),median=median(y_dep),d1=quantile(y_dep,0.1),
                   d9=quantile(y_dep,0.9),var=var(y_dep),var_eps=var(y_eps),var_pred=var(y_pred))]    
  rt$type=name
  rt
}

#' Weighted correclation
#' @export
wt.cor <- function(x,y,w) {
  m1 = sum(x*w)/sum(w)
  v1 = sum((x-m1)^2*w)/sum(w)
  m2 = sum(y*w)/sum(w)
  v2 = sum((y-m2)^2*w)/sum(w)
  cc = sum( (y-m2)*(x-m1)*w)/sum(w)
  return(cc/sqrt(v1*v2))
}

#' Weighted covaraince
#' @export
wt.cov <- function(x,y,w) {
  w = w/sum(w)
  m1 = sum(x*w)
  v1 = sum((x-m1)^2*w)
  m2 = sum(y*w)
  v2 = sum((y-m2)^2*w)
  cc = sum( (y-m2)*(x-m1)*w)
  return(cc)
}

#' Sparse colSums
#' @export
sColSums <- function(M) {
  return(as.numeric((rep(1,dim(M)[1]) %*% M)))
}

#' Sparse rowSums
#' @export
sRowSums <- function(M) {
  return(as.numeric(M %*% rep(1,dim(M)[2])))
}
