# simple simulation according to shimer smith
require(SDMTools)

#' @export
initp <- function(...) {
  p = list()
  p$ay = 0.2
  p$nx = 50
  p$ny = 50
  p$nz = 50
  p$sz = 1
  p$dt  = 4
  p$sep = 0.02
  p$r   = 0.05
  p$rho = -2
  p$Vbar   = 2
  p$b   = 0
  p$c   = 0
  p$alpha = 0.5 # bargaining power
  p$fsize = 50
  p$m   = 0.4
  p$nu  = 0.3 # efficiency of search on the job
  p$K   = 0
  p$lb0  = 0.4
  p$lb1  = 0.1

  #p$pf = function(x,y,z=0,p) x*y + z
  p$pf = function(x,y,z=0,p)  ( x^p$rho +  y^p$rho )^(1/p$rho) + p$ay # haggerdorn law manovski

  sp = list(...)

  for (n in names(sp)) {
    p[[n]] = sp[[n]]
  }

  return(p)
}

#' @export
shimersmith.solve.otj <- function(p,maxiter=400,enditer=maxiter,rs=0.8,rw=0.8,rh=0.5,tol=1e-9) {

  X = spread( (0:(p$nx-1))/p$nx , 2, p$ny )
  Y = spread( (0:(p$ny-1))/p$ny , 1, p$nx )
  #X = exp(qnorm(X))
  #Y = exp(qnorm(Y))

  FF = p$pf(X,Y,0,p)
  W0 = rep(0,p$nx)
  P0 = rep(0,p$ny)
  H  = array(0.8/(p$nx*p$ny),c(p$nx,p$ny))
  U  = rep(0.2/p$nx,p$nx)
  V  = (p$Vbar-0.8) * rep(1/p$ny,p$ny)
  S  = FF / p$r

  t.S    <- tensorFunction(  R[x,y]   ~ b * FF[x,y] - a * W[x] - a * P[y] )
  t.S2   <- tensorFunction(  R[x,y]   ~ (  I( S[x,y2] > S[x,y] ) + 0.5 *  I( S[x,y2] == S[x,y] )  ) * ( a* S[x,y2] - S[x,y] ) *V[y2] )
  t.W0   <- tensorFunction(  R[x]     ~ max(S[x,y],0) * V[y] )
  t.P0   <- tensorFunction(  R[y]     ~ max(S[x,y],0) * U[x] )
  t.P02  <- tensorFunction(  R[y]     ~ ( I(S[x,y]>S[x,y2]) + 0.5*I(S[x,y]==S[x,y2])  ) * S[x,y]* H[x,y] )
  t.H    <- tensorFunction(  R[x,y]   ~ I( S[x,y] >=0 ) * U[x] * V[y] )
  t.H.j2j.in  <- tensorFunction(  R[x,y]   ~ ( I(S[x,y]>S[x,y2]) + 0.5*I(S[x,y]==S[x,y2])  ) * H[x,y2])
  t.H.j2j.out <- tensorFunction(  R[x,y]   ~ ( I(S[x,y]<S[x,y2]) + 0.5*I(S[x,y]==S[x,y2])  ) * V[y2] )
  t.wage <- tensorFunction(   R[x,y] ~ max(0, S[x,y2] - S[x,y] ) * V[y2])

  msg =" all good"
  #plot(S[5,])

  for (i in 1:maxiter) {

    # we compute the intermediary distributions
    H2  = (1 - p$sep) * H # here we might have people with negative surplus, apply S>0
    #U2  = spread(1/p$nx,1,p$nx)   - rowSums(H2)
    #V2  = spread(p$Vbar/p$ny,1,p$ny) - colSums(H2)
    U2 = U + p$sep * rowSums(H)
    V2 = V + p$sep * colSums(H)

    U2 = U2/sum(U2)
    V2 = V2/sum(V2)
    H2 = H2/sum(H2)

    # update mu0 and mu1
    mu0 = p$lb0 *  (p$sep + (1-p$sep)*sum(U)) / (p$sep + (1-p$sep)*sum(V))
    mu1 = p$lb1 *  (1-p$sep) * (1-sum(U) )    / (p$sep + (1-p$sep)*sum(V))
    if (mu1>1) { msg = " mu1 is bigger than 1, increase Vbar!"; break}

    S2 = (t.S(FF,W0,P0,1+p$r,p$r*(1-p$sep)) + p$lb1 * (1 - p$sep) * t.S2(S,V2,p$alpha)  + (1-p$sep) * S)/(1+p$r)
    dS = mean( (S - S2 )^2 )/(.0001 + mean( (S)^2 ))
    if (!is.finite(dS) ) { msg = "NAN in S2"; break;   }
    if (i<enditer) S = rs * S + (1-rs) * S2;
    if (any(!is.finite(S2))) { msg = "NAN in S2"; break;   }
    #lines(S[5,])

    W02 = (  p$b + p$lb0 * p$alpha * t.W0(S,V2) +  W0 ) / (1 + p$r)
    if (i<enditer) W0  = rw * W0 + (1-rw) * W02;
    P02 = (  - p$c  +  mu0 * (1-p$alpha) * t.P0(S,U2)  +  mu1 * (1-p$alpha) * t.P02(S,H2)   + P0   ) / ( 1+ p$r)
    if (i<enditer) P0  = rw * P0 + (1-rw) * P02;
    if (any(is.nan(W0))) { msg = "NAN in W0"; break;   }
    if (any(is.nan(P0))) { msg = "NAN in P0"; break;   }

    # ---------  distributions -------------
    # compute distribution after separations
    # update with meetings
    # we have to remember that worker can search within the period
    # so we compute taking that into consideration
    #u2e     = p$lb0 * (p$sep + (1-p$sep)*sum(U)) * t.H(S,U2,V2)
    #j2j.in  = p$lb1 *  (1-p$sep) * (1-sum(U) ) * t.H.j2j.in(S,H2)  * spread(c(V2),1,p$nx)
    #j2j.out = p$lb1 * t.H.j2j.out(S,V2) * (1 - p$sep) * H
    u2e     = p$lb0 * (p$sep + (1-p$sep)*sum(U)) * t.H(S,U2,V2)
    j2j.in  = p$lb1 *  (1-p$sep) * (1-sum(U) ) * t.H.j2j.in(S,H2)  * spread(c(V2),1,p$nx)
    j2j.out = p$lb1 * t.H.j2j.out(S,V2) * (1 - p$sep) * H
    Hbis    = (1 - p$sep) * H + u2e + j2j.in - j2j.out
    dH = mean( (H - Hbis )^2 )/mean( (H)^2 )

    # updates
    U = U + rowSums(p$sep*H) - rowSums(u2e)
    V = V + colSums(p$sep*H) - colSums(u2e) + colSums(j2j.out) - colSums(j2j.in)

    stopifnot(max( abs(U + rowSums(Hbis) - 1/p$nx))<1e-12)
    stopifnot(max( abs(V + colSums(Hbis) - p$Vbar/p$ny))<1e-12)

    # update distributions using flows
    H  = Hbis #(rh * H + (1-rh)*Hbis)
    # U  = spread(1/p$nx,1,p$nx)   - rowSums(H)
    # V  = spread(p$Vbar/p$ny,1,p$ny) - colSums(H)

    #if ((i %% 100)==0) {
      flog.info("[%3i] dS=%3.3e dH=%3.3e mu0=%3.3e mu1=%3.3e mass=%f",i,dS,dH,mu0,mu1,sum(H) + sum(U))
      #plot(P0)
    #}
    #lines(P0)

    if ( (dS < tol) & (dH < 1e-7)) break;
  }
  flog.info(msg)


  # solving for the wage
  w = p$alpha * (p$sep + p$r) * S + (1-p$sep) * p$r * spread(W0,2,p$ny) -
        p$alpha * (1-p$sep)* p$lb1 * t.wage(S,V2)
  w = w/(1+p$r)

  return(list(S=S,W0=W0,wage=w,H=H,F=FF,V=V,U=U,mu0=mu0,mu1=mu1,W0=W0,P0=P0,V2=V2))
}


#' we simulate a set of movers from shimer smith
#' similar to the jdata
#' @export
shimersmith.simulate.otj <- function(model,p,ni,rsq=0) {

  J1 = array(0,ni)
  J2 = array(0,ni)
  Y1 = array(0,ni)
  Y2 = array(0,ni)
  F1 = array(0,ni)
  F2 = array(0,ni)
  M = array(0,ni)

  # draw worker type in period 1 from the distribution of employed
  X  = sample.int(p$nx,ni,replace = T, pr = rowSums(model$H) )
  fcount = ceiling(ni/(p$ny * p$fsize))

  wsd = sqrt(wt.var( log(c(model$wage)), c(model$H)))
  esd = wsd * sqrt( (1-rsq)/rsq)

  for (i in 1:ni) {
    x = X[i]
    j1 = sample.int(p$ny,1,prob=model$H[x,])
    w1 = log(model$wage[x,j1]) + rnorm(1)*esd
    mover = 0
    j2=-1

    # --------- simulate period 2 ---------

    # draw exogenous separation
    if (rbinom(1,1,prob=p$sep)) {

      mover = -1

      # draw firm from v_{1/2}(y), check matching decision
      if (rbinom(1,1,prob= p$lb0 )) {
        j2 = sample.int(p$ny,1, prob=model$V2)
        # check matching decision
        if (model$S[x,j2]>=0) {
          w2 = log(model$wage[x,j2]) + rnorm(1)*esd
          mover = 2
        } else {
          j2 = 0
          w2 = 0
        }
      } else {
        j2=0
        w2=0
      }

    # no separation  -- search on the job
    } else {
      j2 = j1
      w2 = log(model$wage[x,j2]) + rnorm(1)*esd

      # probability of matching with other firm
      if (rbinom(1,1,prob= p$lb1 )) {
        jt = sample.int(p$ny,1, prob= model$V2)
        # check matching decision
        if (model$S[x,jt] >= model$S[x,j1]) {
          mover = 1
          w2 = log(model$wage[x,jt]) + rnorm(1)*esd
          j2 = jt
        }
      }
    }

    J1[i] = j1
    J2[i] = j2
    Y1[i] = w1
    Y2[i] = w2
    M[i]  = mover
    F1[i] = sample.int(fcount,1,replace = TRUE) + fcount * ( J1[i]-1 )
    if (mover>0) {
      #draw firm, make sure it is a different one
      F2[i] = sample.int(fcount,1,replace = TRUE,prob= !((1:fcount)==F1[i]) ) + fcount * ( J2[i]-1 )

    } else {
      F2[i] = F1[i]
    }
  }

  sdata = data.table(y1=Y1,y2=Y2,t1=J1,t2=J2,f1=F1,f2=F2,x=X,mover=M)

  return(sdata)
}

#' @export
shimersmith.simulateAndEstimate <- function(p,emtrue=FALSE,est_rep=50,est_nbest=1,maxiter=1000) {

  set.seed(218937621)
  flog.info("solving the model")
  model <- shimersmith.solve.otj(p,maxiter=10000,enditer=8000)
  model$wage = model$wage +1  # shift everything up to have well defined logs

  # simulate
  flog.info("simulating data")
  cdata <- shimersmith.simulate.otj(model,p,5e5,rsq=0.9)

  # we cluster the data
  flog.info("cluster data")
  cdata[,f1:=paste(f1)][,f2:=paste(f2)]
  setnames(cdata,"x","k")
  cdata$x=1
  sim = list(sdata=cdata[mover==0],jdata=cdata[mover>0])

  # compute stats on the model
  cstats = sim$sdata[, c(as.list(quantile(y1,seq(0.1,0.9,l=9)) ), list(m=mean(y1),sd=sd(y1))),list(j1=t1)]

  ms    = grouping.getMeasures(sim,"ecdf",Nw=20,y_var = "y1")
  grps  = grouping.classify.once(ms,k = 10,nstart = 1000,iter.max = 200,step=100)
  sim   = grouping.append(sim,grps$best_cluster)

  # we estimate mini-model
  flog.info("linear projection using true types")
  vdec_true = lin.proj(cdata,y_col="y1",k_col="k",j_col = "t1",do.unc = FALSE)
  flog.info("estimating the mini-model")
  res_mini = m2.mini.estimate(sim$jdata,sim$sdata,method = "prof")

  # we estimate
  flog.info("estimating the mixture model")
  ctrl      = em.control(nplot=50,tol=1e-7,dprior=1.001,fixb=TRUE,
                         sd_floor=1e-7,posterior_reg=1e-8,
                         est_rep=est_rep,est_nbest=est_nbest,sdata_subsample=0.1,ncat=25,maxiter=maxiter)
  res_mixt = m2.mixt.estimate.all(sim,nk=6,ctrl)

  #  ------------ counter-factuals ------------
  rme  = sample.stats(sim$sdata,"y1","t1")
  # randomly allocate workers to firms in the cross-section
  sdata.nosorting = copy(sim$sdata)
  sdata.nosorting[,k2 := k[ rank(runif(.N))]]
  sdata.nosorting[,y1_imp := log(model$wage[k2,t1]),list(k2,t1)]
  rme = rbind(rme,sample.stats(sdata.nosorting,"y1_imp","t1"))

  sdata.nosorting = m2.mixt.impute.stayers(sim$sdata,res_mixt$model)
  rme = rbind(rme,sample.stats(sdata.nosorting,"y1_imp","j1"))
  sdata.nosorting[,k2     := k_imp[ rank(runif(.N))]]
  sdata.nosorting[,y1_imp2 := res_mixt$model$A1[j1,k2],list(k2,j1)]
  rme = rbind(rme,sample.stats(sdata.nosorting,"y1_imp2","j1"))

  return(list(model=model,p=p,vdec_true=vdec_true,res_mixt=res_mixt,res_mini=res_mini,mean_effect = rme ,cstats=cstats))
}


#' @export
shimersmith.plot.em <- function(res,model,cdata) {

  # reorder the data
  wm_hat = res$wm
  ID_hat = cdata[,mean(y1),j1][order(j1)][,order(V1)]
  wm_hat = wm_hat[,ID_hat]

  # reorder the model
  wm_true = log(model$wage)
  ID_true = cdata[,mean(y1),t1][order(t1)][,order(V1)]
  wm_true = wm_true[,ID_true]

  ggplot(melt(wm_hat,c('k','l')),aes(x=l,y=value,group=k)) + geom_line(aes(color=factor(k))) + theme_bw() + geom_line(data = melt(wm_true,c('k','l')), linetype=2)
}



