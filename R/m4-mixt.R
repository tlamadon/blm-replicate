# here we want to check a bunch of properties for the EM steps
# model1 and model2 should be 2 consecutive steps
m4.mixt.check <- function(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,...) {

  change = list(...)

  # compute posterior for model1
  res1 = with(model1,{
    taum = array(0,c(Nm,nk))
    lpm  = array(0,c(Nm,nk))
    likm = 0
    for (i in 1:Nm) {
      ltau = log(pk1[JJm[i],])
      lnorm2 = lognormpdf(Y2m[i], A2ma[J1m[i],] + A2mb[J2m[i]]      , S2m[J1m[i],])
      lnorm1 = lognormpdf(Y1m[i] - B12*(Y2m[i] - A2ma[J1m[i],])    , A12[J1m[i],], S12[J1m[i],])
      lnorm3 = lognormpdf(Y3m[i] - B32m*(Y2m[i] - A2ma[J1m[i],] - A2mb[J2m[i]]), A3ma[J2m[i],] + A3mb[J1m[i]], S3m[J2m[i],])
      lnorm4 = lognormpdf(Y4m[i] - B43 *(Y3m[i] - A3ma[J2m[i],])               , A43[J2m[i],], S43[J2m[i],])
      lall = ltau + lnorm2 + lnorm4 + lnorm1  + lnorm3
      lpm[i,] = lall
      if (any(is.na(lall))) { catf("lall is na (%f,%f,%f,%f,%f)!\n",ltau,lnorm1,lnorm2,lnorm3,lnorm4); stop =T; break;}
      likm = likm + logsumexp(lall)
      taum[i,] = exp(lall - logsumexp(lall))
    }

    # compute prior
    lik_prior = (dprior-1) * sum(log(pk1)) # dirichlet distribution
    lik = likm + lik_prior
    list(taum = taum, lpm =lpm, lik=likm,lik_prior=lik_prior,post=lik)
  })

  model2 = copy(model1)
  model2[names(change)] = change[names(change)]

  # compute posterior for model2
  res2 = with(model2,{
    taum = array(0,c(Nm,nk))
    lpm  = array(0,c(Nm,nk))
    likm = 0
    for (i in 1:Nm) {
      ltau = log(pk1[JJm[i],])
      lnorm2 = lognormpdf(Y2m[i], A2ma[J1m[i],] + A2mb[J2m[i]]      , S2m[J1m[i],])
      lnorm1 = lognormpdf(Y1m[i] - B12*(Y2m[i] - A2ma[J1m[i],])    , A12[J1m[i],], S12[J1m[i],])
      lnorm3 = lognormpdf(Y3m[i] - B32m*(Y2m[i] - A2ma[J1m[i],] - A2mb[J2m[i]]), A3ma[J2m[i],] + A3mb[J1m[i]], S3m[J2m[i],])
      lnorm4 = lognormpdf(Y4m[i] - B43 *(Y3m[i] - A3ma[J2m[i],])               , A43[J2m[i],], S43[J2m[i],])
      lall = ltau + lnorm2 + lnorm4 + lnorm1  + lnorm3
      lpm[i,] = lall
      if (any(is.na(lall))) { catf("lall is na (%f,%f,%f,%f,%f)!\n",ltau,lnorm1,lnorm2,lnorm3,lnorm4); stop =T; break;}
      likm = likm + logsumexp(lall)
      taum[i,] = exp(lall - logsumexp(lall))
    }

    # compute prior
    lik_prior = (dprior-1) * sum(log(pk1)) # dirichlet distribution
    lik = likm + lik_prior
    list(taum = taum, lpm =lpm, lik=likm,lik_prior=lik_prior,post=lik)
  })

  # aboid actual 0
  res1$taum[res1$taum==0] = min(res1$taum[res1$taum>0])
  res2$taum[res2$taum==0] = min(res2$taum[res2$taum>0])
  # do the analysis, Evaluate Q(theta | theta^t) , Q(theta^t | theta^t), H(theta | theta^t) and H(theta^t | theta^t)
  Q1 = sum( ( (res1$taum) * res1$lpm ))
  Q2 = sum( ( (res1$taum) * res2$lpm ))
  H1 = - sum( (res1$taum) * log(res1$taum))
  H2 = - sum( (res1$taum) * log(res2$taum))

  if (is.nan(H1)) browser();

  warn_str=""
  test = TRUE
  if (( Q2<Q1) | (H2<H1)) {
    warn_str = "!!!!!!!!!";
    test=FALSE
  }
  flog.info("[emcheck] %s Qd=%4.4e Hd=%4.4e %s",paste(names(change),collapse = ","),  Q2-Q1,H2-H1,warn_str)
  return(test)
}

# here we want to check a bunch of properties for the EM steps
# model1 and model2 should be 2 consecutive steps
m4.mixt.check.stayers <- function(Y1s,Y2s,Y3s,Y4s,J1s,nk,Ns,model1,...) {

  change = list(...)

  # compute posterior for model1
  res1 = with(model1,{
    taum = array(0,c(Ns,nk))
    lpm  = array(0,c(Ns,nk))
    likm = 0
    for (i in 1:Ns) {
      ltau = log(pk0[J1s[i],])
      lnorm2 = lognormpdf(Y2s[i]                                 , A2s[J1s[i],], S2s[J1s[i],])
      lnorm3 = lognormpdf(Y3s[i] - B32s*(Y2s[i] - A2s[J1s[i],])  , A3s[J1s[i],], S3s[J1s[i],])
      lnorm1 = lognormpdf(Y1s[i] - B12 *(Y2s[i] - A2ma[J1s[i],]) , A12[J1s[i],], S12[J1s[i],])
      lnorm4 = lognormpdf(Y4s[i] - B43 *(Y3s[i] - A3ma[J1s[i],]) , A43[J1s[i],], S43[J1s[i],])
      lall = ltau + lnorm2 + lnorm4 + lnorm1  + lnorm3
      lpm[i,] = lall
      if (any(is.na(lall))) { catf("lall is na (%f,%f,%f,%f,%f)!\n",ltau,lnorm1,lnorm2,lnorm3,lnorm4); stop =T; break;}
      likm = likm + logsumexp(lall)
      taum[i,] = exp(lall - logsumexp(lall))
    }

    # compute prior
    lik_prior = (dprior-1) * sum(log(pk1)) # dirichlet distribution
    lik = likm + lik_prior
    list(taum = taum, lpm =lpm, lik=likm,lik_prior=lik_prior,post=lik)
  })

  model2 = copy(model1)
  model2[names(change)] = change[names(change)]

  # compute posterior for model2
  res2 = with(model2,{
    taum = array(0,c(Ns,nk))
    lpm  = array(0,c(Ns,nk))
    likm = 0
    for (i in 1:Ns) {
      ltau = log(pk0[J1s[i],])
      lnorm2 = lognormpdf(Y2s[i]                                 , A2s[J1s[i],], S2s[J1s[i],])
      lnorm3 = lognormpdf(Y3s[i] - B32s*(Y2s[i] - A2s[J1s[i],])  , A3s[J1s[i],], S3s[J1s[i],])
      lnorm1 = lognormpdf(Y1s[i] - B12 *(Y2s[i] - A2ma[J1s[i],]) , A12[J1s[i],], S12[J1s[i],])
      lnorm4 = lognormpdf(Y4s[i] - B43 *(Y3s[i] - A3ma[J1s[i],]) , A43[J1s[i],], S43[J1s[i],])
      lall = ltau + lnorm2 + lnorm4 + lnorm1  + lnorm3
      lpm[i,] = lall
      if (any(is.na(lall))) { catf("lall is na (%f,%f,%f,%f,%f)!\n",ltau,lnorm1,lnorm2,lnorm3,lnorm4); stop =T; break;}
      likm = likm + logsumexp(lall)
      taum[i,] = exp(lall - logsumexp(lall))
    }

    # compute prior
    lik_prior = (dprior-1) * sum(log(pk1)) # dirichlet distribution
    lik = likm + lik_prior
    list(taum = taum, lpm =lpm, lik=likm,lik_prior=lik_prior,post=lik)
  })

  # do the analysis, Evaluate Q(theta | theta^t) , Q(theta^t | theta^t), H(theta | theta^t) and H(theta^t | theta^t)
  Q1 = sum( ( (res1$taum) * res1$lpm ))
  Q2 = sum( ( (res1$taum) * res2$lpm ))
  H1 = - sum( (res1$taum) * log(res1$taum))
  H2 = - sum( (res1$taum) * log(res2$taum))

  warn_str=""
  test = TRUE
  if (( Q2<Q1) | (H2<H1)) {
    warn_str = "!!!!!!!!!";
    test=FALSE
  }
  catf("[emcheck] %s Qd=%4.4e Hd=%4.4e %s\n",paste(names(change),collapse = ","),  Q2-Q1,H2-H1,warn_str)
  return(test)
}

# here we want to check a bunch of properties for the EM steps
# model1 and model2 should be 2 consecutive steps
m4.mixt.check.stayers.unc <- function(Y1s,Y2s,Y3s,Y4s,J1s,nk,Ns,model1,...) {

  change = list(...)

  # compute posterior for model1
  res1 = with(model1,{
    taum = array(0,c(Ns,nk))
    lpm  = array(0,c(Ns,nk))
    likm = 0
    for (i in 1:Ns) {
      ltau = log(pk0[J1s[i],])
      lnorm2 = lognormpdf(Y2s[i]                                 , A2s[J1s[i],], S2s[J1s[i],])
      lnorm3 = lognormpdf(Y3s[i] - B32s*(Y2s[i] - A2s[J1s[i],])  , A3s[J1s[i],], S3s[J1s[i],])
      lnorm1 = lognormpdf(Y1s[i] - B12 *(Y2s[i] - A2s[J1s[i],])  , A12[J1s[i],], S12[J1s[i],])
      lnorm4 = lognormpdf(Y4s[i] - B43 *(Y3s[i] - A3s[J1s[i],])  , A43[J1s[i],], S43[J1s[i],])
      lall = ltau + lnorm2 + lnorm4 + lnorm1  + lnorm3
      lpm[i,] = lall
      if (any(is.na(lall))) { catf("lall is na (%f,%f,%f,%f,%f)!\n",ltau,lnorm1,lnorm2,lnorm3,lnorm4); stop =T; break;}
      likm = likm + logsumexp(lall)
      taum[i,] = exp(lall - logsumexp(lall))
    }

    # compute prior
    lik_prior = (dprior-1) * sum(log(pk1)) # dirichlet distribution
    lik = likm + lik_prior
    list(taum = taum, lpm =lpm, lik=likm,lik_prior=lik_prior,post=lik)
  })

  model2 = copy(model1)
  model2[names(change)] = change[names(change)]

  # compute posterior for model2
  res2 = with(model2,{
    taum = array(0,c(Ns,nk))
    lpm  = array(0,c(Ns,nk))
    likm = 0
    for (i in 1:Ns) {
      ltau = log(pk0[J1s[i],])
      lnorm2 = lognormpdf(Y2s[i]                                 , A2s[J1s[i],], S2s[J1s[i],])
      lnorm3 = lognormpdf(Y3s[i] - B32s*(Y2s[i] - A2s[J1s[i],])  , A3s[J1s[i],], S3s[J1s[i],])
      lnorm1 = lognormpdf(Y1s[i] - B12 *(Y2s[i] - A2s[J1s[i],])  , A12[J1s[i],], S12[J1s[i],])
      lnorm4 = lognormpdf(Y4s[i] - B43 *(Y3s[i] - A3s[J1s[i],])  , A43[J1s[i],], S43[J1s[i],])
      lall = ltau + lnorm2 + lnorm4 + lnorm1  + lnorm3
      lpm[i,] = lall
      if (any(is.na(lall))) { catf("lall is na (%f,%f,%f,%f,%f)!\n",ltau,lnorm1,lnorm2,lnorm3,lnorm4); stop =T; break;}
      likm = likm + logsumexp(lall)
      taum[i,] = exp(lall - logsumexp(lall))
    }

    # compute prior
    lik_prior = (dprior-1) * sum(log(pk1)) # dirichlet distribution
    lik = likm + lik_prior
    list(taum = taum, lpm =lpm, lik=likm,lik_prior=lik_prior,post=lik)
  })

  # do the analysis, Evaluate Q(theta | theta^t) , Q(theta^t | theta^t), H(theta | theta^t) and H(theta^t | theta^t)
  Q1 = sum( ( (res1$taum) * res1$lpm ))
  Q2 = sum( ( (res1$taum) * res2$lpm ))
  H1 = - sum( (res1$taum) * log(res1$taum))
  H2 = - sum( (res1$taum) * log(res2$taum))

  warn_str=""
  test = TRUE
  if (( Q2<Q1) | (H2<H1)) {
    warn_str = "!!!!!!!!!";
    test=FALSE
  }
  catf("[emcheck] %s Qd=%4.4e Hd=%4.4e %s\n",paste(names(change),collapse = ","),  Q2-Q1,H2-H1,warn_str)
  return(test)
}


m4.mixt.new.from.mini <- function(nk,model0) {
  NNm   = model0$Nm
  NNs   = model0$Ns
  jdata = model.mini4.simulate.movers(model0,NNm)
  sdata = model.mini4.simulate.stayers(model0,NNs)
  nf =model0$nf
  model = m4.mixt.new(nk,nf)
  jdata[,k := ceiling(rank(alpha)/.N*nk)]

  ctrl   = em.control(nplot=20,model0=model)
  n=jdata[,.N]
  tau = array(0.01,c(n,nk))
  tau[1:n + n*(jdata$k-1)]=1-(nk-1)*0.01
  model$B12 = model0$r1
  model$B43 = model0$r4
  model$B32m = model0$rho32m
  res    = m4.mixt.rhoext.movers(  jdata, sdata, model,   em.control(ctrl,maxiter=1,check_lik=F,fixb=T,tau=tau),em.order = 1)
  return(res$model)
}

# create a random model for EM with
# endogenous mobility with multinomial pr
#' @export
m4.mixt.new <-function(nk,nf,inc=F,from.mini=NA) {

  model = list()
  # model for Y1|Y2,l,k for movers and stayes
  model$A12    = array(0.9*(1 + 0.5*rnorm(nf*nk)),c(nf,nk))
  model$B12    = 0.4 + 0.5*runif(1)
  model$S12    = array(1+0.5*runif(nf*nk),c(nf,nk))
  # model for Y4|Y3,l,k for movers and stayes
  model$A43    = array(0.9*(1 + 0.5*rnorm(nf*nk)),c(nf,nk))
  model$B43    = 0.4 + 0.5*runif(1)
  model$S43    = array(1+0.5*runif(nf*nk),c(nf,nk))
  # model for p(K | l ,l') for movers
  model$pk1    = rdirichlet(nf*nf,rep(1,nk))
  # model for p(K | l ) for stayers
  model$pk0    = rdirichlet(nf,rep(1,nk))

  # model for Y2 | k,l,l'   Y2 = N( \mu_{kl} + \xi_{l'} , \sigma_{kl} )
  # for movers
  model$A2ma      = array(0.9*(1 + 0.5*rnorm(nf*nk)),c(nf,nk))
  model$A2mb      = c(0,rnorm( nf-1 ))
  model$S2m       = array(1+0.5*runif(nf*nk)  , c( nf, nk  ))
  # model for Y3 | Y2, k,l,l'   Y2 = N( \mu_{kl} + \xi_{l'} , \sigma_{kl} )
  # for movers
  model$A3ma      = array(0.9*(1 + 0.5*rnorm(nf*nk)),c(nf,nk))
  model$A3mb      = c(0,rnorm( nf -1 ))
  model$S3m       = array(1+0.5*runif(nf*nk)  , c( nf, nk  ))
  model$B32m      = 0.4 + 0.3*runif(1)

  # model for Y2 | k,l,l'   Y2 = N( \mu_{kl} + \xi_{l'} , \sigma_{kl} )
  # for movers
  model$A2s       = array(0.9*(1 + 0.5*rnorm(nf*nk)),c(nf,nk))
  for (l in 1:nf) {model$A2s[l,] = sort(model$A2s[l,])}
  model$S2s       = array(1+0.5*runif(nf*nk)  , c( nf, nk  ))
  # model for Y3 | Y2, k,l,l'   Y2 = N( \mu_{kl} + \xi_{l'} , \sigma_{kl} )
  # for movers
  model$A3s       = array(0.9*(1 + 0.5*rnorm(nf*nk)),c(nf,nk))
  model$S3s       = array(1+0.5*runif(nf*nk)  , c( nf, nk  ))
  model$B32s      = 0.4 + 0.5*runif(1)

  if (inc) {
    for (l in 1:nf) {
      model$A12[l,] = sort(model$A12[l,])
      model$A43[l,] = sort(model$A43[l,])
      model$A2ma[l,] = sort(model$A2ma[l,])
      model$A3ma[l,] = sort(model$A3ma[l,])
    }
  }

  model$nk    = nk
  model$nf    = nf
  model$likm  =-Inf
  return(model)
}

# ---- SIMULATING ----


#' Simulates data (movers and stayers) and attached firms ids. Firms have all same expected size.
#' @export
m4.mixt.simulate.sim <- function(model,fsize,smult=1,mmult=1) {
  jdata       = m4.mixt.simulate.movers(model,model$NNm*mmult)
  sdata       = m4.mixt.simulate.stayers(model,model$NNs*smult)

  # create some firm ids
  sdata <- sdata[,f1:=paste("F",j1 + model$nf*(sample.int(.N/fsize,.N,replace=T)-1),sep=""),j1]
  sdata <- sdata[,j1b:=j1]
  sdata <- sdata[,j1true := j1]
  jdata <- jdata[,j1true := j1][,j2true := j2]
  jdata <- jdata[,j1c:=j1]
  jdata <- jdata[,f1:=sample( unique(sdata[j1b %in% j1c,f1]) ,.N,replace=T),j1c]
  jdata <- jdata[,j2c:=j2]
  jdata <- jdata[,f2:=sample( unique(sdata[j1b %in% j2c,f1])  ,.N,replace=T),j2c]
  jdata$j2c=NULL
  jdata$j1c=NULL
  sdata$j1b=NULL
  sdata[,f2:=f1]

  sim = list(sdata=sdata,jdata=jdata)
  return(sim)
}


#' @export
m4.mixt.simulate.movers <- function(model,NNm) {

  J1 = array(0,sum(NNm))
  J2 = array(0,sum(NNm))
  Y0 = array(0,sum(NNm))
  Y1 = array(0,sum(NNm))
  Y2 = array(0,sum(NNm))
  Y3 = array(0,sum(NNm))
  Y4 = array(0,sum(NNm))
  Y5 = array(0,sum(NNm))
  K  = array(0,sum(NNm))

  A12  = model$A12
  B12  = model$B12
  S12  = model$S12
  A43  = model$A43
  B43  = model$B43
  S43  = model$S43
  pk1  = model$pk1
  A2ma = model$A2ma
  A2mb = model$A2mb
  S2m  = model$S2m
  A3ma = model$A3ma
  A3mb = model$A3mb
  S3m  = model$S3m
  B32m = model$B32m
  nk   = model$nk
  nf   = model$nf

  i =1
  for (l1 in 1:nf) for (l2 in 1:nf) {
    I = i:(i+NNm[l1,l2]-1)
    ni = length(I)
    jj = l1 + nf*(l2 -1)
    J1[I] = l1
    J2[I] = l2

    # draw k
    Ki = sample.int(nk,ni,T,pk1[jj,])
    K[I] = Ki

    # draw Y2, Y3
    Y2[I]  = A2ma[l1,Ki] + A2mb[l2] + S2m[l1,Ki] * rnorm(ni)
    Y3[I]  = A3ma[l2,Ki] + A3mb[l1] + S3m[l2,Ki] * rnorm(ni) + B32m*(Y2[I] - A2ma[l1,Ki] - A2mb[l2]) # here we include A2mb

    # draw Y1, Y4
    Y1[I] = rnorm(ni)*S12[l1,Ki] + A12[l1,Ki] + B12*(Y2[I] - A2ma[l1,Ki]) # the equation for Y1|Y2, no A2mb here
    Y4[I] = rnorm(ni)*S43[l2,Ki] + A43[l2,Ki] + B43*(Y3[I] - A3ma[l2,Ki]) # the equation for Y4|Y3, no A3mb here

    i = i + NNm[l1,l2]
  }

  jdatae = data.table(k=K,y1=Y1,y2=Y2,y3=Y3,y4=Y4,j1=J1,j2=J2)
  return(jdatae)
}

#' @export
m4.mixt.simulate.stayers <- function(model,NNs) {

  J1 = array(0,sum(NNs))
  Y1 = array(0,sum(NNs))
  Y2 = array(0,sum(NNs))
  Y3 = array(0,sum(NNs))
  Y4 = array(0,sum(NNs))
  K  = array(0,sum(NNs))

  A12  = model$A12
  B12  = model$B12
  S12  = model$S12
  A43  = model$A43
  B43  = model$B43
  S43  = model$S43
  pk0  = model$pk0
  A2s  = model$A2s
  A2ma = model$A2ma
  S2s  = model$S2s
  A3s  = model$A3s
  S3s  = model$S3s
  A3ma = model$A3ma
  B32s = model$B32s
  nk   = model$nk
  nf   = model$nf

  i =1
  for (l1 in 1:nf) {
    I = i:(i+NNs[l1]-1)
    ni = length(I)
    J1[I] = l1

    # draw k
    Ki = sample.int(nk,ni,T,pk0[l1,])
    K[I] = Ki

    # draw Y2, Y3
    Y2[I]  = A2s[l1,Ki] + S2s[l1,Ki] * rnorm(ni)
    Y3[I]  = A3s[l1,Ki] + S3s[l1,Ki] * rnorm(ni) + B32s*(Y2[I] - A2s[l1,Ki])

    # draw Y1, Y4
    Y1[I] = rnorm(ni)*S12[l1,Ki] + A12[l1,Ki] + B12*(Y2[I] - A2ma[l1,Ki]) # the equation for Y1|Y2, same as movers!
    Y4[I] = rnorm(ni)*S43[l1,Ki] + A43[l1,Ki] + B43*(Y3[I] - A3ma[l1,Ki]) # the equation for Y4|Y3, same as movers!

    i = i + NNs[l1]
  }

  sdatae = data.table(k=K,y1=Y1,y2=Y2,y3=Y3,y4=Y4,j1=J1)
  return(sdatae)
}

#' @export
m4.mixt.impute.movers <- function(jdatae,model) {

  A12  = model$A12
  B12  = model$B12
  S12  = model$S12
  A43  = model$A43
  B43  = model$B43
  S43  = model$S43
  pk1  = model$pk1
  A2ma = model$A2ma
  A2mb = model$A2mb
  S2m  = model$S2m
  A3ma = model$A3ma
  A3mb = model$A3mb
  S3m  = model$S3m
  B32m = model$B32m
  nk   = model$nk
  nf   = model$nf

  # ------  impute K, Y1, Y4 on jdata ------- #
  jdatae.sim = copy(jdatae)
  jdatae.sim[, c('k_imp','y1_imp','y2_imp','y3_imp','y4_imp') := {
    ni = .N
    jj = j1 + nf*(j2-1)
    Ki  = sample.int(nk,.N,prob = pk1[jj,],replace=T)
    # draw Y2, Y3
    Y2  = A2ma[j1,Ki] + A2mb[j2] + S2m[j1,Ki] * rnorm(ni)
    Y3  = A3ma[j2,Ki] + A3mb[j1] + S3m[j2,Ki] * rnorm(ni) + B32m*(Y2 - A2ma[j1,Ki] - A2mb[j2])
    # draw Y1, Y4
    Y1 = rnorm(ni)*S12[j1,Ki] + A12[j1,Ki] + B12*(Y2 - A2ma[j1,Ki]) # the equation for Y1|Y2, no A2mb here
    Y4 = rnorm(ni)*S43[j2,Ki] + A43[j2,Ki] + B43*(Y3 - A3ma[j2,Ki]) # the equation for Y1|Y2, no A2mb here
    list(Ki,Y1,Y2,Y3,Y4)
  },list(j1,j2)]

  return(jdatae.sim)
}


#' @export
m4.mixt.impute.stayers <- function(sdatae,model) {

  A12  = model$A12
  B12  = model$B12
  S12  = model$S12
  A43  = model$A43
  B43  = model$B43
  S43  = model$S43
  pk0  = model$pk0
  A2s  = model$A2s
  S2s  = model$S2s
  A3s  = model$A3s
  S3s  = model$S3s
  B32s = model$B32s
  A2ma = model$A2ma
  A3ma = model$A3ma
  nk   = model$nk
  nf   = model$nf

  # ------  impute K, Y1, Y4 on jdata ------- #
  sdatae.sim = copy(sdatae)
  sdatae.sim[, c('k_imp','y1_imp','y2_imp','y3_imp','y4_imp') := {
    ni = .N
    Ki  = sample.int(nk,.N,prob = pk0[j1,],replace=T)
    # draw Y2, Y3
    Y2  = A2s[j1,Ki] + S2s[j1,Ki] * rnorm(ni)
    Y3  = A3s[j1,Ki] + S3s[j1,Ki] * rnorm(ni) + B32s*(Y2 - A2s[j1,Ki])
    # draw Y1, Y4
    Y1 = rnorm(ni)*S12[j1,Ki] + A12[j1,Ki] + B12*(Y2 - A2ma[j1,Ki]) # the equation for Y1|Y2, same as movers!
    Y4 = rnorm(ni)*S43[j1,Ki] + A43[j1,Ki] + B43*(Y3 - A3ma[j1,Ki]) # the equation for Y1|Y2, same as movers!
    list(Ki,Y1,Y2,Y3,Y4)
  },list(j1)]

  return(sdatae.sim)
}

m4.mixt.simulatebest <- function() {

  # load the grid
  load("../figures/src/em-endo-full_rhogrid-halton-6x10.dat",verbose=F)

  # find the best
  dd = data.frame()
  for (r in rr) {
    dd = rbind(dd,data.frame(rho=r$model$B32m,lik=r$lik,time=r$time.taken,step=r$step,dlik=r$dlik))
  }
  rbest = rr[[which.max(dd$lik)]]
  cat(sprintf("%i evaluations, best value is %f\n",length(rr),rbest$lik))

  # get number of movers
  load("../figures/src/em-endo-info.dat",verbose=F)

  # reweight the statyers to 30,0000
  tot = NNs[,sum(ni)]
  NNs[,ni := round(ni * 30000 /sum(ni)) ]
  setkey(NNs,j)
  NNs = NNs[,ni]

  # get the movers matrix
  NNm = acast(NNm,j1~j2,value.var="ni")

  # ----- simulate ------ #
  nk = rbest$model$nk;
  nf = rbest$model$nf;
  model = rbest$model
  jdatae = m4.mixt.full.simulate.movers(model,NNm)
  sdatae = m4.mixt.full.simulate.stayers(model,NNs)
  jdatae[,m:=1][,w:=1]
  sdatae[,m:=0][,w:=tot/.N]
  sdatae[,j2:=j1]
  adatae = rbind(sdatae,sdatae)
  cat(sprintf("simulated data with %i stayers and %i movers \n",sdatae[,.N],jdatae[,.N]))
  return(adatae)
}


#' @export
makePosteriorStochastic <- function(tau,m) {
  if (m==0) return(tau)
  # we compute the frequency from M draw from each
  N    = dim(tau)[1]
  nk   = dim(tau)[2]
  tau2 = array(0,c(N,m))
  tau3 = array(0,c(N,nk))
  for (i in 1:N) {
    tau2[i,] = sample.int(nk,m,replace=T,prob = tau[i,])
  }
  for (k in 1:nk) {
    tau3[,k] = rowSums(tau2==k)
  }
  tau3 = tau3/rowSums(tau3)
  return(tau3)
}

# ------- EM ESTIMATION -----

#' Estimates both movers and stayers with several starting values.
#'
#' 1) the rhos are estaimted using the stayers
#' 2) the parameters for movers are estimated using movers only:
#'   a) starts with forced parallel values
#'   b) spread the parallel values, keep them fix, estimate proportions
#'   c) relax parallel assumption, but do not estimate effect of past firms
#'   d) remove all restrictions
#' 3) repeat step 2 est_rep times, select est_nbest best likelihoods
#' 4) compute worst connectedness among types, pick starting values with best worst value
#' 5) estimate parameters for movers:
#'   a) use subamble
#'   b) force AKM restriction and estimate
#'   c) relax AKM restriction and estimate
#'
#' The function can make use of the cl parameter (a created cluster) to perform
#' the estimation in parallel across nodes.
#'
#' @param cl a cluster created with makeCluster from the parallel package
#'
#' @export
m4.mixt.estimate.all <-  function(sim,nk=6,ctrl,cl=NA) {

  start.time <- Sys.time()
  nf = max(sim$sdata$j1);

  # get overall mean and sd of y2 (this is for exploration in different starting values)
  mm = mean(sim$sdata$y2)
  ms = 2*sd(sim$sdata$y2)

  # extracting the rho values using the stayers
  sim$sdata[,y1:=y1_bu]
  res_rhos             = m4.mini.getvar.stayers.unc2.opt(sim$sdata,diff=ctrl$rho_in_diff)

  # prepare the control function  for the EM estimation
  ctrl    = em.control(ctrl=ctrl,est_rho=c(FALSE,TRUE,FALSE,FALSE,TRUE,FALSE))

  nf = max(sim$sdata$j1);
  model_start      = m4.mixt.new(nk,nf)
  model_start$B12  = res_rhos$r1
  model_start$B43  = res_rhos$r4
  model_start$B32m = 0.6

  # evaluate parallel estimates
  res_para = m4.mixt.rhoext.movers(sim$jdata, sim$sdata,model_start,ctrl=em.control(ctrl,cstr_type="para",textapp="para0",fixb=FALSE))

  # use cluster if available
  if (!any(is.na(cl))) {
    flog.info("cluster -- exporting objects to nodes")
    # export environment to nodes
    clusterExport(cl,c("res_para","sim","ctrl"),environment())
    mylapply <- function(...) parLapply(cl,...)
    nnodes=length(cl)
  } else {
    mylapply <- function(...) lapply(...)
    nnodes=1
  }

  flog.info("starting repetitions with %i nodes",nnodes)
  rr  = mylapply(1:ctrl$est_rep, function(i) {
    res_mixt = list()
    tryCatch({
      res_para$model$A2mb[]=0
      res_para$model$A3mb[]=0
      res_para$model$A2ma  = spread(sort(rnorm(nk))*ms+mm,1,nf)
      res_para$model$A3ma  = res_para$model$A2ma
      res_para$model$A12   = res_para$model$A2ma
      res_para$model$A43   = res_para$model$A2ma

      # -------- unconstrained ------- #
      res_para_fixm = m4.mixt.rhoext.movers(sim$jdata,sim$sdata,res_para$model,ctrl=em.control(ctrl,fixb=FALSE,cstr_type="para", textapp=sprintf("paraf (%i/%i)",i,ctrl$est_rep),fixm=T))
      res_para_new  = m4.mixt.rhoext.movers(sim$jdata,sim$sdata,res_para_fixm$model,ctrl=em.control(ctrl,fixb=FALSE,textapp=sprintf("para1 (%i/%i)",i,ctrl$est_rep),cstr_type="para"))
      res_mixt_noAmb= m4.mixt.rhoext.movers(sim$jdata,sim$sdata,res_para_new$model,ctrl=em.control(ctrl,textapp=sprintf("move1 (%i/%i)",i,ctrl$est_rep),fixb=T,est_Amb=F))
      res_mixt      = m4.mixt.rhoext.movers(sim$jdata,sim$sdata,res_mixt_noAmb$model,ctrl=em.control(ctrl,textapp=sprintf("move2 (%i/%i)",i,ctrl$est_rep),fixb=T,est_Amb=T))

      # ------ compute connectedness ----- #
      res_mixt$connectedness = model.connectiveness(res_mixt$model)
      res_mixt$rep_id = i
    }, error = function(e) {flog.error("error in rep %i!",i);print(e)})
    flog.info("done with reptitions %i/%i",i,ctrl$est_rep)
    res_mixt
  })

  # extract likelihoods and connectedness
  rrd = ldply(rr,function(r) {
    data.frame(lik_mixt = r$model$likm,connectedness = r$connectedness,i=r$rep_id)
  })

  # selecting best starting value
  rrd     = data.table(rrd)
  rrd[, sel:=-1]
  rrd.sub = rrd[order(-lik_mixt)][1:ctrl$est_nbest]
  rrd[i %in% rrd.sub$i, sel:=0]
  Ibest = rrd.sub[order(-connectedness)][1,i]
  res_mixt = rr[[Ibest]]
  rrd[i==Ibest, sel:=1]

  # sub-sample the stayers for computational reasons (if too large)
  if (ctrl$sdata_subredraw==TRUE) {
    sim$sdata[,sample := rank(runif(.N))/.N<=ctrl$sdata_subsample,j1]
  }

  # we run the stayers
  sim$sdata[,y1:=y1_bu]
  res_mixt = m4.mixt.rhoext.stayers(sim$jdata, sim$sdata[sample==1], res_mixt$model, em.control(ctrl,cstr_type="akm",cstr_val=0.05,textapp="stay1",fixb=F))
  sim$sdata[,y1:=y1_bu]
  res_mixt = m4.mixt.rhoext.stayers(sim$jdata, sim$sdata[sample==1], res_mixt$model, em.control(ctrl,cstr_type="none",textapp="stay2"))
  res_mixt$second_stage_reps     = rrd
  res_mixt$second_stage_reps_all = rr

  # ------ compute linear decomposition ------- #
  stayer_share = sum(res_mixt$model$NNs)/(sum(res_mixt$model$NNm)*ctrl$sdata_subsample + sum(res_mixt$model$NNs))
  proj_unc  = m4.mixt.vdec(res_mixt$model,ctrl$vdec_sim_size,stayer_share,"y2")
  res_mixt$vdec = proj_unc

  end.time <- Sys.time()
  res_mixt$time.taken <- end.time - start.time

  return(res_mixt)
}


#' Computes the variance decomposition by simulation
#' @export
m4.mixt.vdec <- function(model,nsim,stayer_share=1,ydep="y2") {

  # simulate movers/stayers, and combine
  NNm = model$NNm
  NNs = model$NNs
  NNm[!is.finite(NNm)]=0
  NNs[!is.finite(NNs)]=0

  NNs = round(NNs*nsim*stayer_share/sum(NNs))
  NNm = round(NNm*nsim*(1-stayer_share)/sum(NNm))
  flog.info("computing var decomposition with ns=%i nm=%i",sum(NNs),sum(NNm))


  # we simulate from the model both movers and stayers
  sdata.sim = m4.mixt.simulate.stayers(model,NNs)
  jdata.sim = m4.mixt.simulate.movers(model,NNm)
  sdata.sim = rbind(sdata.sim[,list(j1,k,y2)],jdata.sim[,list(j1,k,y2)])
  proj_unc  = lin.proj(sdata.sim,ydep,"k","j1");

  return(proj_unc)
}

#' Compute mean effects
#' @export
m4.mixt.meaneffect <- function(model) {
  NNs = model$NNs*10 # used 10% sample
  NNm = model$NNm
  share_s  = sum(NNs)/(sum(NNm) + sum(NNs))
  share_m  = sum(NNm)/(sum(NNm) + sum(NNs))

  NNs = round(NNs*1e6*share_s/sum(NNs))
  NNm = round(NNm*1e6*share_m/sum(NNm))

  # we simulate from the model both movers and stayers
  sdata = m4.mixt.simulate.stayers(model,NNs)
  jdata = m4.mixt.simulate.movers(model,NNm)
  sdata = rbind(sdata[,list(j1,k,y2)],jdata[,list(j1,k,y2)])

  # compute decomposition
  #vdec  = lin.proj(sdata,y_col = "y1",k_col="k",j_col = "j1")
  #res_bs$mixt_all[[nn]]$vdec_1m = vdec
  rt    = sample.stats(sdata,"y2","j1","pk0")

  # then we set the distribution to uniform
  model_nosort     = copy(model)
  model_nosort$pk0 = spread(colSums(model_nosort$pk0 * spread(NNs/(sum(NNs)),2,model_nosort$nk)),1,model_nosort$nf)

  # movers
  dpk1 = m2.get.pk1(model)
  pk   = dpk1[,pr_k[1],k][,V1]
  model_nosort$pk1      = spread(pk,1,model$nf * model$nf)

  # simulate from uniform
  sdata = m4.mixt.simulate.stayers(model_nosort,NNs)
  jdata = m4.mixt.simulate.movers(model_nosort,NNm)
  sdata = rbind(sdata[,list(j1,k,y2)],jdata[,list(j1,k,y2)])
  rt2   = sample.stats(sdata,"y2","j1","pku")

  return(rbind(rt,rt2))
}

#' Estimates the dynamic model parameters, but keeps the rhos fixed
#'
#' @export
m4.mixt.rhoext.movers <- function(jdatae,sdatae,model,ctrl,em.order=1) {

  start.time <- Sys.time()
  tic <- tic.new()

  dprior = ctrl$dprior
  model0 = ctrl$model0
  taum = ctrl$tau

  ### ----- GET MODEL  ---
  nk   = model$nk
  nf   = model$nf
  # Y1 | Y2
  A12  = model$A12
  B12  = model$B12
  S12  = model$S12
  A43  = model$A43
  B43  = model$B43
  S43  = model$S43
  ## --movers
  pk1  = model$pk1
  A2ma = model$A2ma
  A2mb = model$A2mb
  S2m  = model$S2m
  A3ma = model$A3ma
  A3mb = model$A3mb
  A2ma = model$A2ma
  A2mb = model$A2mb
  S3m  = model$S3m
  B32m = model$B32m

  # ----- GET DATA
  # movers
  Y1m = jdatae$y1
  Y2m = jdatae$y2
  Y3m = jdatae$y3
  Y4m = jdatae$y4
  J1m = jdatae$j1
  J2m = jdatae$j2
  JJm = J1m + nf*(J2m-1)
  Nm = jdatae[,.N]

  # get the constraints
  CS12 = cons.pad(cons.get(ctrl$cstr_type[1],ctrl$cstr_val[1],nk,nf),nk*nf*0, nk*nf*3 + (nf-1)*2)
  CS2  = cons.pad(cons.get(ctrl$cstr_type[2],ctrl$cstr_val[2],nk,nf),nk*nf*1, nk*nf*2 + (nf-1)*2)
  CS32 = cons.pad(cons.get(ctrl$cstr_type[3],ctrl$cstr_val[3],nk,nf),nk*nf*2, nk*nf*1 + (nf-1)*2)
  CS43 = cons.pad(cons.get(ctrl$cstr_type[4],ctrl$cstr_val[4],nk,nf),nk*nf*3, nk*nf*0 + (nf-1)*2)
  # combine them
  CS = cons.bind(cons.bind(cons.bind(CS12,CS2),CS32),CS43)

  # create the stationary contraints
  if (ctrl$fixb==T) {
    CS2 = cons.pad(cons.fixb(nk,nf,4),0,(nf-1)*2)
    CS  = cons.bind(CS2,CS)
  }

  # create a constraint for the variances
  if (ctrl$model_var==T) {
    CSw = cons.none(nk,nf*4)
  } else{
    CS12 = cons.pad(cons.mono_k(nk,nf),nk*nf*0, nk*nf*3)
    CS2  = cons.pad(cons.mono_k(nk,nf),nk*nf*1, nk*nf*2)
    CS32 = cons.pad(cons.mono_k(nk,nf),nk*nf*2, nk*nf*1)
    CS43 = cons.pad(cons.mono_k(nk,nf),nk*nf*3, nk*nf*0)
    CSw  = cons.bind(cons.bind(cons.bind(CS12,CS2),CS32),CS43)
    CSw$meq = length(CSw$H)
  }

  # prepare matrices aggregated at the type level
  Dkj1f     = diag(nf)  %x% rep(1,nf) %x% diag(nk)            # A[k,l] coefficients for j1
  Dkj2f     = rep(1,nf) %x% diag(nf)  %x% diag(nk)            # A[k,l] coefficients for j2
  Dj1f_bis  = (diag(nf)  %x% rep(1,nf) %x% rep(1,nk))[,2:nf]  # C[l]   coefficients for j1
  Dj2f_bis  = (rep(1,nf) %x% diag(nf)  %x% rep(1,nk))[,2:nf]  # C[l]   coefficients for j2

  # regression matrix for the variance
  XXV = rBind(
      cBind(    Dkj1f, 0*Dkj1f,  0*Dkj2f, 0*Dkj2f ),
      cBind(  0*Dkj1f,   Dkj1f,  0*Dkj2f, 0*Dkj2f ),
      cBind(  0*Dkj1f, 0*Dkj1f,    Dkj2f, 0*Dkj2f ),
      cBind(  0*Dkj1f, 0*Dkj1f,  0*Dkj2f,   Dkj2f )
  )

  ## --- prepare regressions covariates --- #
  # create the depend variables

  lik_old  = -Inf
  lik      = -Inf
  lik_best = -Inf
  liks = 0
  likm=0

  lpt1 = array(0,c(Nm,nk))
  lpt2 = array(0,c(Nm,nk))
  lpt3 = array(0,c(Nm,nk))
  lpt4 = array(0,c(Nm,nk))
  lp   = array(0,c(Nm,nk))


  tic("prep")
  stop = F;
  for (step in 1:ctrl$maxiter) {

    model1 = list(nk=nk,nf=nk,A12=A12,B12=B12,S12=S12,
                  A2ma=A2ma,A2mb=A2mb,S2m=S2m,
                  A3ma=A3ma,A3mb=A3mb,S3m=S3m,B32m=B32m,
                  A43=A43,B43=B43,S43=S43,
                  pk1=pk1,dprior=dprior)

    ### ---------- E STEP ------------- #
    # compute the tau probabilities and the likelihood
    if (is.na(taum[1]) | (step>1)) {

      # for efficiency we want to group by (l1,l2)
      for (l1 in 1:nf) for (l2 in 1:nf) {
        I = which( (J1m==l1) & (J2m==l2))
        ll = l1 + nf*(l2-1)
        if (length(I)==0) next;

        for (k in 1:nk) {
          lpt2[I,k] = lognormpdf(Y2m[I]                                        , A2ma[l1,k] + A2mb[l2] , S2m[l1,k])
          lpt1[I,k] = lognormpdf(Y1m[I] - B12 *(Y2m[I] - A2ma[l1,k])           , A12[l1,k], S12[l1,k])
          lpt3[I,k] = lognormpdf(Y3m[I] - B32m*(Y2m[I] - A2ma[l1,k] - A2mb[l2]), A3ma[l2,k] + A3mb[l1], S3m[l2,k])
          lpt4[I,k] = lognormpdf(Y4m[I] - B43 *(Y3m[I] - A3ma[l2,k])           , A43[l2,k], S43[l2,k])

          # sum the log of the periods
          lp[I,k] = log(pk1[ll,k]) + lpt1[I,k] + lpt2[I,k] + lpt3[I,k] + lpt4[I,k]
        }
      }

      liks     = sum(logRowSumExp(lp))
      taum     = exp(lp - spread(logRowSumExp(lp),2,nk)) # normalize the k probabilities Pr(k|Y1,Y2,Y3,Y4,l)

      # compute prior
      lik_prior = (dprior-1) * sum(log(pk1))
      lik = liks + lik_prior

    } else {
      cat("skiping first max step, using supplied posterior probabilities\n")
    }

    tic("estep")

    if (stop) break;

    # ---------- MAX STEP ------------- #
    # taum = makePosteriorStochastic(tau = taum,m = ctrl$stochastic) # if we want to implement stochastic EM

    # we start by recovering the posterior weight, and the variances for each term
    rwm  = c(t(taum + ctrl$posterior_reg))
    wv12 = c(t(S12[J1m,]))^2
    wv2  = c(t(S2m[J1m,]))^2
    wv32 = c(t(S3m[J2m,]))^2
    wv43 = c(t(S43[J2m,]))^2

    # here we switch between updating the rho paramerters
    # and updating the other "level" parameters
    update_levels = ((em.order + (step))%%2 )==0
    if (all(ctrl$est_rho[1:3]==FALSE)) update_levels=TRUE;

    if (update_levels==TRUE) {
      DYY  = array(0,c(nk,nf,nf,4))
      WWT  = array(1e-7,c(nk,nf,nf,4))

      for (l1 in 1:nf) for (l2 in 1:nf) {
        I = which( (J1m==l1) & (J2m==l2))
        if (length(I)==0) next;
        for (k in 1:nk) {
          # compute the posterior weight, it's not time specific
          ww = sum(taum[I,k]+ ctrl$posterior_reg)

          # construct dependent for each time period k,l2,l1,
          DYY[k,l2,l1,1] = sum(  (Y1m[I] - B12  * Y2m[I]) * (taum[I,k]+ ctrl$posterior_reg) )/ww
          DYY[k,l2,l1,2] = sum(  (Y2m[I]                ) * (taum[I,k]+ ctrl$posterior_reg) )/ww
          DYY[k,l2,l1,3] = sum(  (Y3m[I] - B32m * Y2m[I]) * (taum[I,k]+ ctrl$posterior_reg) )/ww
          DYY[k,l2,l1,4] = sum(  (Y4m[I] - B43  * Y3m[I]) * (taum[I,k]+ ctrl$posterior_reg) )/ww

          # Scaling the weight by the time specific variance
          WWT[k,l2,l1,1] = ww/pmax(ctrl$sd_floor,S12[l1,k]^2)
          WWT[k,l2,l1,2] = ww/pmax(ctrl$sd_floor,S2m[l1,k]^2)
          WWT[k,l2,l1,3] = ww/pmax(ctrl$sd_floor,S3m[l2,k]^2)
          WWT[k,l2,l1,4] = ww/pmax(ctrl$sd_floor,S43[l2,k]^2)
        }
      }

     # modifiy the regressor to take into account the possibly changed B12, B23m, B43
     XXf     = rBind(
      cBind(    Dkj1f,  -B12*Dkj1f,       0*Dkj1f,      0*Dkj1f,      0*Dj2f_bis,  0*Dj1f_bis),
      cBind(  0*Dkj1f,       Dkj1f,       0*Dkj1f,      0*Dkj1f,        Dj2f_bis,  0*Dj1f_bis),
      cBind(  0*Dkj1f, -B32m*Dkj1f,         Dkj2f,      0*Dkj1f,  -B32m*Dj2f_bis,    Dj1f_bis),
      cBind(  0*Dkj1f,     0*Dkj1f,    -B43*Dkj2f,        Dkj2f,      0*Dj2f_bis,  0*Dj1f_bis)
     )

     if (ctrl$fixm==TRUE) {
       fit = c(as.numeric(t(A12)),as.numeric(t(A2ma)),as.numeric(t(A3ma)),as.numeric(t(A43)),A2mb[2:nf],A3mb[2:nf])
     } else {
       fit   = slm.wfitc(XXf,as.numeric(DYY),as.numeric(WWT),CS,scaling=0)$solution
       A12   = t(rdim(fit[1:(nk*nf)],nk,nf))
       A2ma  = t(rdim(fit[(nk*nf+1):(2*nk*nf)],nk,nf))
       A3ma  = t(rdim(fit[(2*nk*nf+1):(3*nk*nf)],nk,nf))
       A43   = t(rdim(fit[(3*nk*nf+1):(4*nk*nf)],nk,nf))
       A2mb  = c(0,fit[(4*nk*nf+1):(4*nk*nf+nf-1)])
       A3mb  = c(0,fit[(4*nk*nf+nf):(4*nk*nf+2*(nf-1))])
       if (ctrl$check_lik) {
         m4.mixt.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,A12=A12,A2ma=A2ma,A3ma=A3ma,A43=A43,A2mb=A2mb,A3mb=A3mb)
       }
     }

    # compute the variances!!!!
    DYY_bar   = array(0,c(nk,nf,nf,4))
    DYY_bar[] = XXf%*%fit
    DYYV      = array(0,c(nk,nf,nf,4))

    for (l1 in 1:nf) for (l2 in 1:nf) {
      I = which( (J1m==l1) & (J2m==l2))
      if (length(I)==0) next;
      for (k in 1:nk) {
        # construct dependent for each time period k,l2,l1,
        ww = sum(taum[I,k]+ ctrl$posterior_reg)
        DYYV[k,l2,l1,1] = sum(  (Y1m[I] - B12  * Y2m[I] - DYY_bar[k,l2,l1,1])^2 * (taum[I,k]+ ctrl$posterior_reg) )/ww
        DYYV[k,l2,l1,2] = sum(  (Y2m[I]                 - DYY_bar[k,l2,l1,2])^2 * (taum[I,k]+ ctrl$posterior_reg) )/ww
        DYYV[k,l2,l1,3] = sum(  (Y3m[I] - B32m * Y2m[I] - DYY_bar[k,l2,l1,3])^2 * (taum[I,k]+ ctrl$posterior_reg) )/ww
        DYYV[k,l2,l1,4] = sum(  (Y4m[I] - B43  * Y3m[I] - DYY_bar[k,l2,l1,4])^2 * (taum[I,k]+ ctrl$posterior_reg) )/ww
      }
    }

     fitv  = slm.wfitc(XXV,as.numeric(DYYV),as.numeric(WWT),CSw)$solution
     S12   = sqrt(t(rdim(fitv[1:(nk*nf)],nk,nf)))
     S2m   = sqrt(t(rdim(fitv[(nk*nf+1):(2*nk*nf)],nk,nf)))
     S3m   = sqrt(t(rdim(fitv[(2*nk*nf+1):(3*nk*nf)],nk,nf)))
     S43   = sqrt(t(rdim(fitv[(3*nk*nf+1):(4*nk*nf)],nk,nf)))
     if (ctrl$check_lik) {
        m4.mixt.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,A12=A12,A2ma=A2ma,A3ma=A3ma,A43=A43,A2mb=A2mb,A3mb=A3mb,S12=S12,S2m=S2m,S3m=S3m,S43=S43)
     }
     S12[S12<ctrl$sd_floor]=ctrl$sd_floor # having a variance of exacvtly 0 creates problem in the likelihood
     S2m[S2m<ctrl$sd_floor]=ctrl$sd_floor
     S3m[S3m<ctrl$sd_floor]=ctrl$sd_floor # having a variance of exacvtly 0 creates problem in the likelihood
     S43[S43<ctrl$sd_floor]=ctrl$sd_floor

    } else {
      # estimate the rhos -- use the other parameters
      if (ctrl$est_rho[1]==TRUE) {
        fit = lm.wfit( as.matrix(c(t(Y2m - A2ma[J1m,]))), c(t(Y1m - A12[J1m,])),rwm/wv12)
        B12 = coefficients(fit)
        if (ctrl$check_lik) { m4.mixt.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,B12=B12)}
      }
      if (ctrl$est_rho[3]==TRUE) {
        fit = lm.wfit( as.matrix(c(t(Y3m - A3ma[J2m,]))), c(t(Y4m - A43[J2m,])),rwm/wv43)
        B43 = coefficients(fit)
        if (ctrl$check_lik) { m4.mixt.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,B43=B43)}
      }
      if (ctrl$est_rho[2]==TRUE) {
        fit = lm.wfit( as.matrix(c(t(Y2m - A2ma[J1m,] - A2mb[J2m]))), c(t(Y3m - A3ma[J2m,] - A3mb[J1m])),rwm/wv32)
        B32m = coefficients(fit)
        if (ctrl$check_lik) { m4.mixt.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,B32m=B32m)}
      }
    }

    tic("mstep-ols")

    ## -------- PK probabilities ------------ #
    ## --- movers --- #
    for (l1 in 1:nf) for (l2 in 1:nf) {
      jj = l1 + nf*(l2-1)
      I = which(JJm == jj)
      if (length(I)>1) {
        pk1[jj,] = colSums(taum[I,])
      } else if (length(I)==0) { # this deals with the case where the cell is empty
        pk1[jj,] = 1/nk
      } else {
        pk1[jj,] = taum[I,]
      }
      pk1[jj,] = (pk1[jj,] + dprior-1 )/(sum(pk1[jj,] + dprior -1 ))
    }
    if (ctrl$check_lik) { m4.mixt.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,pk1=pk1)}

    #check_lik = computeLik(Y1m,Y2m,Y3m,Y4m,A12,B12,S12,A43,B43,S43,A2ma,A2mb,S2m,A3ma,A3mb,B32m,S3m)
    #if (check_lik<lik) cat("lik did not go down on pk1 update\n")

    tic("mstep-pks")

    # checking model fit
    if ((!any(is.na(model0))) & ((step %% ctrl$nplot) == (ctrl$nplot-1))) {
      I1 = order(colSums(A2ma))
      I2 = order(colSums(model0$A2ma))
      rr = addmom(A2ma[,I1],model0$A2ma[,I2],"A2ma")
      rr = addmom(A2mb,model0$A2mb,"A2mb",rr)
      rr = addmom(A3ma[,I1],model0$A3ma[,I2],"A3ma",rr)
      rr = addmom(A3mb,model0$A3mb,"A3mb",rr)
      rr = addmom(A12[,I1],model0$A12[,I2],"A12",rr)
      rr = addmom(A43[,I1],model0$A43[,I2],"A43",rr)
      rr = addmom(c(B12,B32m,B43), c(model0$B12,model0$B32m,model0$B43), "rhos", rr,type="rho")
      rr = addmom(S2m[,I1], model0$S2m[,I2], "S2m", rr,type="var")
      rr = addmom(S3m[,I1], model0$S3m[,I2], "S3m", rr,type="var")
      rr = addmom(S12[,I1],model0$S12[,I2],"S12",rr,type="var")
      rr = addmom(S43[,I1],model0$S43[,I2],"S43",rr,type="var")

      rr = addmom(pk1,model0$pk1,"pk1",rr,type="pr")

      print(ggplot(rr,aes(x=val2,y=val1,color=type)) + geom_point() + facet_wrap(~name,scale="free") + theme_bw() + geom_abline(linetype=2))
    } else {
      if ((step %% ctrl$nplot) == (ctrl$nplot-1)) {
        mm = list(A12=A12,A2ma=A2ma,A3ma=A3ma,A43=A43)
        m4.mixt.movers.plotw(mm)
      }
    }

    # -------- check convergence ------- #
    dlik = (lik - lik_old)/abs(lik_old)
    lik_old = lik
    lik_best = pmax(lik_best,lik)
    if ( (step %% ctrl$ncat) == 0) flog.info("[%3i][%s] levels=%i lik=%4.4f dlik=%4.4e r12=%.3f r23=%.3f r43=%.3f",step,ctrl$textapp,update_levels+0,lik,dlik,B12,B32m,B43);
    if (step>10) if (abs(dlik)<ctrl$tol) break;

    tic("loop-wrap")
  }

  flog.info("[%3i][%s][final] levels=%i lik=%4.4f dlik=%4.4e r12=%.3f r23=%.3f r43=%.3f",step,ctrl$textapp,update_levels+0,lik,dlik,B12,B32m,B43);

  # Y1 | Y2
  model$A12  = A12
  model$B12  = B12
  model$S12  = S12
  model$A43  = A43
  model$B43  = B43
  model$S43  = S43
  ## --movers --
  model$pk1  = pk1
  model$A2ma = A2ma
  model$A2mb = A2mb
  model$S2m  = S2m
  model$A3ma = A3ma
  model$A3mb = A3mb
  model$S3m  = S3m
  model$B32m = B32m

  model$NNm = acast(jdatae[,.N,list(j1,j2)],j1~j2,fill=0,value.var="N")
  model$likm=liks

  end.time <- Sys.time()
  time.taken <- end.time - start.time

  return(list(tic = tic(), model=model,lik=lik,step=step,dlik=dlik,time.taken=time.taken,ctrl=ctrl,liks=liks,likm=likm))
}

#' This is a version of the EM with a perturbation. This is to
#' explore the likelihood surface.
#'
#' @export
m4.mixt.rhoext.movers.pert <- function(jdatae,sdatae,model,ctrl,em.order=1) {

  start.time <- Sys.time()
  tic <- tic.new()

  dprior = ctrl$dprior
  model0 = ctrl$model0
  taum = ctrl$tau

  ### ----- GET MODEL  ---
  nk   = model$nk
  nf   = model$nf
  # Y1 | Y2
  A12  = model$A12
  B12  = model$B12
  S12  = model$S12
  A43  = model$A43
  B43  = model$B43
  S43  = model$S43
  ## --movers
  pk1  = model$pk1
  A2ma = model$A2ma
  A2mb = model$A2mb
  S2m  = model$S2m
  A3ma = model$A3ma
  A3mb = model$A3mb
  A2ma = model$A2ma
  A2mb = model$A2mb
  S3m  = model$S3m
  B32m = model$B32m

  # ----- GET DATA
  # movers
  Y1m = jdatae$y1
  Y2m = jdatae$y2
  Y3m = jdatae$y3
  Y4m = jdatae$y4
  J1m = jdatae$j1
  J2m = jdatae$j2
  JJm = J1m + nf*(J2m-1)
  Nm = jdatae[,.N]

  # get the constraints
  CS12 = cons.pad(cons.get(ctrl$cstr_type[1],ctrl$cstr_val[1],nk,nf),nk*nf*0, nk*nf*3 + (nf-1)*2)
  CS2  = cons.pad(cons.get(ctrl$cstr_type[2],ctrl$cstr_val[2],nk,nf),nk*nf*1, nk*nf*2 + (nf-1)*2)
  CS32 = cons.pad(cons.get(ctrl$cstr_type[3],ctrl$cstr_val[3],nk,nf),nk*nf*2, nk*nf*1 + (nf-1)*2)
  CS43 = cons.pad(cons.get(ctrl$cstr_type[4],ctrl$cstr_val[4],nk,nf),nk*nf*3, nk*nf*0 + (nf-1)*2)
  # combine them
  CS = cons.bind(cons.bind(cons.bind(CS12,CS2),CS32),CS43)

  # create the stationary contraints
  if (ctrl$fixb==T) {
    CS2 = cons.pad(cons.fixb(nk,nf,4),0,(nf-1)*2)
    CS  = cons.bind(CS2,CS)
  }

  # create a constraint for the variances
  if (ctrl$model_var==T) {
    CSw = cons.none(nk,nf*4)
  } else{
    CS12 = cons.pad(cons.mono_k(nk,nf),nk*nf*0, nk*nf*3)
    CS2  = cons.pad(cons.mono_k(nk,nf),nk*nf*1, nk*nf*2)
    CS32 = cons.pad(cons.mono_k(nk,nf),nk*nf*2, nk*nf*1)
    CS43 = cons.pad(cons.mono_k(nk,nf),nk*nf*3, nk*nf*0)
    CSw  = cons.bind(cons.bind(cons.bind(CS12,CS2),CS32),CS43)
    CSw$meq = length(CSw$H)
  }

  # prepare matrices aggregated at the type level
  Dkj1f     = diag(nf)  %x% rep(1,nf) %x% diag(nk)            # A[k,l] coefficients for j1
  Dkj2f     = rep(1,nf) %x% diag(nf)  %x% diag(nk)            # A[k,l] coefficients for j2
  Dj1f_bis  = (diag(nf)  %x% rep(1,nf) %x% rep(1,nk))[,2:nf]  # C[l]   coefficients for j1
  Dj2f_bis  = (rep(1,nf) %x% diag(nf)  %x% rep(1,nk))[,2:nf]  # C[l]   coefficients for j2

  # regression matrix for the variance
  XXV = rBind(
    cBind(    Dkj1f, 0*Dkj1f,  0*Dkj2f, 0*Dkj2f ),
    cBind(  0*Dkj1f,   Dkj1f,  0*Dkj2f, 0*Dkj2f ),
    cBind(  0*Dkj1f, 0*Dkj1f,    Dkj2f, 0*Dkj2f ),
    cBind(  0*Dkj1f, 0*Dkj1f,  0*Dkj2f,   Dkj2f )
  )

  ## --- prepare regressions covariates --- #
  # create the depend variables

  lik_old  = -Inf
  lik      = -Inf
  lik_best = -Inf
  liks = 0
  likm=0

  lpt1 = array(0,c(Nm,nk))
  lpt2 = array(0,c(Nm,nk))
  lpt3 = array(0,c(Nm,nk))
  lpt4 = array(0,c(Nm,nk))
  lp   = array(0,c(Nm,nk))


  tic("prep")
  stop = F;
  for (step in 1:ctrl$maxiter) {

    model1 = list(nk=nk,nf=nk,A12=A12,B12=B12,S12=S12,
                  A2ma=A2ma,A2mb=A2mb,S2m=S2m,
                  A3ma=A3ma,A3mb=A3mb,S3m=S3m,B32m=B32m,
                  A43=A43,B43=B43,S43=S43,
                  pk1=pk1,dprior=dprior)

    ### ---------- E STEP ------------- #
    # compute the tau probabilities and the likelihood
    if (is.na(taum[1]) | (step>1)) {

      # for efficiency we want to group by (l1,l2)
      for (l1 in 1:nf) for (l2 in 1:nf) {
        I = which( (J1m==l1) & (J2m==l2))
        ll = l1 + nf*(l2-1)
        if (length(I)==0) next;

        for (k in 1:nk) {
          lpt2[I,k] = lognormpdf(Y2m[I]                                        , A2ma[l1,k] + A2mb[l2] , S2m[l1,k])
          lpt1[I,k] = lognormpdf(Y1m[I] - B12 *(Y2m[I] - A2ma[l1,k])           , A12[l1,k], S12[l1,k])
          lpt3[I,k] = lognormpdf(Y3m[I] - B32m*(Y2m[I] - A2ma[l1,k] - A2mb[l2]), A3ma[l2,k] + A3mb[l1], S3m[l2,k])
          lpt4[I,k] = lognormpdf(Y4m[I] - B43 *(Y3m[I] - A3ma[l2,k])           , A43[l2,k], S43[l2,k])

          # sum the log of the periods
          lp[I,k] = log(pk1[ll,k]) + lpt1[I,k] + lpt2[I,k] + lpt3[I,k] + lpt4[I,k]
        }
      }

      liks     = sum(logRowSumExp(lp))
      taum     = exp(lp - spread(logRowSumExp(lp),2,nk)) # normalize the k probabilities Pr(k|Y1,Y2,Y3,Y4,l)

      # compute prior
      lik_prior = (dprior-1) * sum(log(pk1))
      lik = liks + lik_prior

    } else {
      cat("skiping first max step, using supplied posterior probabilities\n")
    }

    tic("estep")

    if (stop) break;

    # ---------- MAX STEP ------------- #
    # taum = makePosteriorStochastic(tau = taum,m = ctrl$stochastic) # if we want to implement stochastic EM

    # we start by recovering the posterior weight, and the variances for each term
    rwm  = c(t(taum + ctrl$posterior_reg))
    wv12 = c(t(S12[J1m,]))^2
    wv2  = c(t(S2m[J1m,]))^2
    wv32 = c(t(S3m[J2m,]))^2
    wv43 = c(t(S43[J2m,]))^2

    # here we switch between updating the rho paramerters
    # and updating the other "level" parameters
    update_levels = ((em.order + (step))%%2 )==0
    if (all(ctrl$est_rho[1:3]==FALSE)) update_levels=TRUE;

    # in this perturbation we randomly weight data points
    WD = rep(1,Nm)
    if (ctrl$pert_type=="data") {
      WD = as.numeric(rdirichlet(1,rep(ctrl$pert_beta,Nm)))*Nm
      taum = taum*spread(WD,2,nk)
    }
    if (ctrl$pert_type=="sa") {
      taum = taum^(ctrl$pert_beta)
      taum = taum/spread(rowSums(taum),2,nk)
    }

    if (update_levels==TRUE) {
      DYY  = array(0,c(nk,nf,nf,4))
      WWT  = array(1e-7,c(nk,nf,nf,4))

      for (l1 in 1:nf) for (l2 in 1:nf) {
        I = which( (J1m==l1) & (J2m==l2))
        if (length(I)==0) next;
        for (k in 1:nk) {
          # compute the posterior weight, it's not time specific
          ww = sum(taum[I,k]+ ctrl$posterior_reg)

          # construct dependent for each time period k,l2,l1,
          DYY[k,l2,l1,1] = sum(  (Y1m[I] - B12  * Y2m[I]) * (taum[I,k]+ ctrl$posterior_reg) )/ww
          DYY[k,l2,l1,2] = sum(  (Y2m[I]                ) * (taum[I,k]+ ctrl$posterior_reg) )/ww
          DYY[k,l2,l1,3] = sum(  (Y3m[I] - B32m * Y2m[I]) * (taum[I,k]+ ctrl$posterior_reg) )/ww
          DYY[k,l2,l1,4] = sum(  (Y4m[I] - B43  * Y3m[I]) * (taum[I,k]+ ctrl$posterior_reg) )/ww

          # Scaling the weight by the time specific variance
          WWT[k,l2,l1,1] = ww/pmax(ctrl$sd_floor,S12[l1,k]^2)
          WWT[k,l2,l1,2] = ww/pmax(ctrl$sd_floor,S2m[l1,k]^2)
          WWT[k,l2,l1,3] = ww/pmax(ctrl$sd_floor,S3m[l2,k]^2)
          WWT[k,l2,l1,4] = ww/pmax(ctrl$sd_floor,S43[l2,k]^2)
        }
      }

      # modifiy the regressor to take into account the possibly changed B12, B23m, B43
      XXf     = rBind(
        cBind(    Dkj1f,  -B12*Dkj1f,       0*Dkj1f,      0*Dkj1f,      0*Dj2f_bis,  0*Dj1f_bis),
        cBind(  0*Dkj1f,       Dkj1f,       0*Dkj1f,      0*Dkj1f,        Dj2f_bis,  0*Dj1f_bis),
        cBind(  0*Dkj1f, -B32m*Dkj1f,         Dkj2f,      0*Dkj1f,  -B32m*Dj2f_bis,    Dj1f_bis),
        cBind(  0*Dkj1f,     0*Dkj1f,    -B43*Dkj2f,        Dkj2f,      0*Dj2f_bis,  0*Dj1f_bis)
      )

      if (ctrl$fixm==TRUE) {
        fit = c(as.numeric(t(A12)),as.numeric(t(A2ma)),as.numeric(t(A3ma)),as.numeric(t(A43)),A2mb[2:nf],A3mb[2:nf])
      } else {
        fit   = slm.wfitc(XXf,as.numeric(DYY),as.numeric(WWT),CS,scaling=0)$solution
        A12   = t(rdim(fit[1:(nk*nf)],nk,nf))
        A2ma  = t(rdim(fit[(nk*nf+1):(2*nk*nf)],nk,nf))
        A3ma  = t(rdim(fit[(2*nk*nf+1):(3*nk*nf)],nk,nf))
        A43   = t(rdim(fit[(3*nk*nf+1):(4*nk*nf)],nk,nf))
        A2mb  = c(0,fit[(4*nk*nf+1):(4*nk*nf+nf-1)])
        A3mb  = c(0,fit[(4*nk*nf+nf):(4*nk*nf+2*(nf-1))])
        if (ctrl$check_lik) {
          m4.mixt.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,A12=A12,A2ma=A2ma,A3ma=A3ma,A43=A43,A2mb=A2mb,A3mb=A3mb)
        }
      }

      # compute the variances!!!!
      DYY_bar   = array(0,c(nk,nf,nf,4))
      DYY_bar[] = XXf%*%fit
      DYYV      = array(0,c(nk,nf,nf,4))

      for (l1 in 1:nf) for (l2 in 1:nf) {
        I = which( (J1m==l1) & (J2m==l2))
        if (length(I)==0) next;
        for (k in 1:nk) {
          # construct dependent for each time period k,l2,l1,
          ww = sum(taum[I,k]+ ctrl$posterior_reg)
          DYYV[k,l2,l1,1] = sum(  (Y1m[I] - B12  * Y2m[I] - DYY_bar[k,l2,l1,1])^2 * (taum[I,k]+ ctrl$posterior_reg) )/ww
          DYYV[k,l2,l1,2] = sum(  (Y2m[I]                 - DYY_bar[k,l2,l1,2])^2 * (taum[I,k]+ ctrl$posterior_reg) )/ww
          DYYV[k,l2,l1,3] = sum(  (Y3m[I] - B32m * Y2m[I] - DYY_bar[k,l2,l1,3])^2 * (taum[I,k]+ ctrl$posterior_reg) )/ww
          DYYV[k,l2,l1,4] = sum(  (Y4m[I] - B43  * Y3m[I] - DYY_bar[k,l2,l1,4])^2 * (taum[I,k]+ ctrl$posterior_reg) )/ww
        }
      }

      fitv  = slm.wfitc(XXV,as.numeric(DYYV),as.numeric(WWT),CSw)$solution
      S12   = sqrt(t(rdim(fitv[1:(nk*nf)],nk,nf)))
      S2m   = sqrt(t(rdim(fitv[(nk*nf+1):(2*nk*nf)],nk,nf)))
      S3m   = sqrt(t(rdim(fitv[(2*nk*nf+1):(3*nk*nf)],nk,nf)))
      S43   = sqrt(t(rdim(fitv[(3*nk*nf+1):(4*nk*nf)],nk,nf)))
      if (ctrl$check_lik) {
        m4.mixt.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,A12=A12,A2ma=A2ma,A3ma=A3ma,A43=A43,A2mb=A2mb,A3mb=A3mb,S12=S12,S2m=S2m,S3m=S3m,S43=S43)
      }
      S12[S12<ctrl$sd_floor]=ctrl$sd_floor # having a variance of exacvtly 0 creates problem in the likelihood
      S2m[S2m<ctrl$sd_floor]=ctrl$sd_floor
      S3m[S3m<ctrl$sd_floor]=ctrl$sd_floor # having a variance of exacvtly 0 creates problem in the likelihood
      S43[S43<ctrl$sd_floor]=ctrl$sd_floor

    } else {
      # estimate the rhos -- use the other parameters
      if (ctrl$est_rho[1]==TRUE) {
        fit = lm.wfit( as.matrix(c(t(Y2m - A2ma[J1m,]))), c(t(Y1m - A12[J1m,])),rwm/wv12)
        B12 = coefficients(fit)
        if (ctrl$check_lik) { m4.mixt.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,B12=B12)}
      }
      if (ctrl$est_rho[3]==TRUE) {
        fit = lm.wfit( as.matrix(c(t(Y3m - A3ma[J2m,]))), c(t(Y4m - A43[J2m,])),rwm/wv43)
        B43 = coefficients(fit)
        if (ctrl$check_lik) { m4.mixt.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,B43=B43)}
      }
      if (ctrl$est_rho[2]==TRUE) {
        fit = lm.wfit( as.matrix(c(t(Y2m - A2ma[J1m,] - A2mb[J2m]))), c(t(Y3m - A3ma[J2m,] - A3mb[J1m])),rwm/wv32)
        B32m = coefficients(fit)
        if (ctrl$check_lik) { m4.mixt.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,B32m=B32m)}
      }
    }

    tic("mstep-ols")

    ## -------- PK probabilities ------------ #
    ## --- movers --- #
    for (l1 in 1:nf) for (l2 in 1:nf) {
      jj = l1 + nf*(l2-1)
      I = which(JJm == jj)
      if (length(I)>1) {
        pk1[jj,] = colSums(taum[I,])
      } else if (length(I)==0) { # this deals with the case where the cell is empty
        pk1[jj,] = 1/nk
      } else {
        pk1[jj,] = taum[I,]
      }
      pk1[jj,] = (pk1[jj,] + dprior-1 )/(sum(pk1[jj,] + dprior -1 ))
    }
    if (ctrl$check_lik) { m4.mixt.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,pk1=pk1)}

    #check_lik = computeLik(Y1m,Y2m,Y3m,Y4m,A12,B12,S12,A43,B43,S43,A2ma,A2mb,S2m,A3ma,A3mb,B32m,S3m)
    #if (check_lik<lik) cat("lik did not go down on pk1 update\n")

    tic("mstep-pks")

    # checking model fit
    if ((!any(is.na(model0))) & ((step %% ctrl$nplot) == (ctrl$nplot-1))) {
      I1 = order(colSums(A2ma))
      I2 = order(colSums(model0$A2ma))
      rr = addmom(A2ma[,I1],model0$A2ma[,I2],"A2ma")
      rr = addmom(A2mb,model0$A2mb,"A2mb",rr)
      rr = addmom(A3ma[,I1],model0$A3ma[,I2],"A3ma",rr)
      rr = addmom(A3mb,model0$A3mb,"A3mb",rr)
      rr = addmom(A12[,I1],model0$A12[,I2],"A12",rr)
      rr = addmom(A43[,I1],model0$A43[,I2],"A43",rr)
      rr = addmom(c(B12,B32m,B43), c(model0$B12,model0$B32m,model0$B43), "rhos", rr,type="rho")
      rr = addmom(S2m[,I1], model0$S2m[,I2], "S2m", rr,type="var")
      rr = addmom(S3m[,I1], model0$S3m[,I2], "S3m", rr,type="var")
      rr = addmom(S12[,I1],model0$S12[,I2],"S12",rr,type="var")
      rr = addmom(S43[,I1],model0$S43[,I2],"S43",rr,type="var")

      rr = addmom(pk1,model0$pk1,"pk1",rr,type="pr")

      print(ggplot(rr,aes(x=val2,y=val1,color=type)) + geom_point() + facet_wrap(~name,scale="free") + theme_bw() + geom_abline(linetype=2))
    } else {
      if ((step %% ctrl$nplot) == (ctrl$nplot-1)) {
        mm = list(A12=A12,A2ma=A2ma,A3ma=A3ma,A43=A43)
        m4.mixt.movers.plotw(mm)
      }
    }

    # -------- check convergence ------- #
    dlik = (lik - lik_old)/abs(lik_old)
    lik_old = lik
    lik_best = pmax(lik_best,lik)
    if ( (step %% ctrl$ncat) == 0) flog.info("[%3i][%s] levels=%i lik=%4.4f dlik=%4.4e top=%4.4e r12=%.3f r23=%.3f r43=%.3f",step,ctrl$textapp,update_levels+0,lik,dlik,lik_best,B12,B32m,B43);
    if (step>10) if (abs(dlik)<ctrl$tol) break;

    tic("loop-wrap")
  }

  flog.info("[%3i][%s][final] levels=%i lik=%4.4f dlik=%4.4e r12=%.3f r23=%.3f r43=%.3f",step,ctrl$textapp,update_levels+0,lik,dlik,B12,B32m,B43);

  # Y1 | Y2
  model$A12  = A12
  model$B12  = B12
  model$S12  = S12
  model$A43  = A43
  model$B43  = B43
  model$S43  = S43
  ## --movers --
  model$pk1  = pk1
  model$A2ma = A2ma
  model$A2mb = A2mb
  model$S2m  = S2m
  model$A3ma = A3ma
  model$A3mb = A3mb
  model$S3m  = S3m
  model$B32m = B32m

  model$NNm = acast(jdatae[,.N,list(j1,j2)],j1~j2,fill=0,value.var="N")
  model$likm=liks

  end.time <- Sys.time()
  time.taken <- end.time - start.time

  return(list(tic = tic(), model=model,lik=lik,step=step,dlik=dlik,time.taken=time.taken,ctrl=ctrl,liks=liks,likm=likm))
}

#' EM for endogenous mobility and modeling of full joint Y1,Y2,Y3,Y4
#' in this version, rho1, rho2 and rhotilda are taken as given
#' B12, B34, B32 are scalars
#' @export
m4.mixt.rhoext.stayers <- function(jdatae,sdatae,model,ctrl) {

  start.time <- Sys.time()
  tic <- tic.new()

  dprior = ctrl$dprior
  model0 = ctrl$model0
  tau = ctrl$tau

  ## ----- GET MODEL ---- #
  nk   = model$nk
  nf   = model$nf
  # Y1 | Y2
  A12  = model$A12
  B12  = model$B12
  S12  = model$S12
  A43  = model$A43
  B43  = model$B43
  S43  = model$S43
  A2ma  = model$A2ma
  A3ma  = model$A3ma
  ## -- stayers --
  pk0  = model$pk0
  A2s  = model$A2s
  S2s  = model$S2s
  A3s  = model$A3s
  S3s  = model$S3s
  B32s = model$B32s

  ## ----- GET DATA ---- #
  # stayers
  Y1s = sdatae$y1
  Y2s = sdatae$y2
  Y3s = sdatae$y3
  Y4s = sdatae$y4
  J1s = sdatae$j1
  Ns  = sdatae[,.N]

  S2s[S2s<ctrl$sd_floor]=ctrl$sd_floor # having a variance of exacvtly 0 creates problem in the likelihood
  S3s[S3s<ctrl$sd_floor]=ctrl$sd_floor


  ## ----- PREPARE ---- #
  ## --- stayers --- #
  Dj1 = array(0,c(Ns,nf))
  Dj1[1:Ns + Ns*(J1s-1)]=1
  # duplicate them nk times (to use posterior weigths)
  Dkj1    = kronecker(Dj1,diag(nk))
  Dkj1    = as.matrix.csr(Dkj1)
  Dkj0    = Dkj1*0

  # get the constraints
  CS2  = cons.pad(cons.get(ctrl$cstr_type[2],ctrl$cstr_val[2],nk,nf),0, nk*nf)
  CS32 = cons.pad(cons.get(ctrl$cstr_type[3],ctrl$cstr_val[3],nk,nf), nk*nf,0)
  # combine them
  CS   = cons.bind(CS2,CS32)

  # add the fixb contraints
  if (ctrl$fixb==TRUE) {
    # interactions must be equal to movers
    CSf   = cons.mono_k(nk,nf)
    CSf$H = c(t(model$A12[,2:nk]-model$A12[,1:(nk-1)]))
    CSf$meq = length(CSf$H)

    # save constraints for both time periods
    CSf = cons.bind(cons.pad(CSf,0,nk*nf),cons.pad(CSf,nk*nf,0))
    CS  = cons.bind(CS,CSf)
  }

  # ---- prepare regressions covariates ---- #

  XXV     = rBind(
      cBind(  Dkj1, Dkj0 ),
      cBind(  Dkj0, Dkj1 )
  )

  lik_old  = -Inf
  lik      = -Inf
  lik_best = -Inf
  likm=0
  dlik_ma = 1
  dlik_smooth = 0

  model1 = copy(model)
  model1$dprior = ctrl$dprior

  tic("prep")

  stop = F;
  for (step in 1:ctrl$maxiter) {


    model1[c('A2s','A3s','S3s','S2s','pk0','B32s')] = list(A2s,A3s,S3s,S2s,pk0,B32s)

    # ---------- E STEP -------------#
    # compute the tau probabilities and the likelihood
    if (is.na(tau[1]) | (step>1)) {

      # --- stayers ---
      taus = array(0,c(Ns,nk))
      liks = 0
      for (i in 1:Ns) {
        ltau = log(pk0[J1s[i],])
        lnorm2 = lognormpdf(Y2s[i]                                 , A2s[J1s[i],], S2s[J1s[i],])
        lnorm3 = lognormpdf(Y3s[i] - B32s*(Y2s[i] - A2s[J1s[i],])  , A3s[J1s[i],], S3s[J1s[i],])
        lnorm1 = lognormpdf(Y1s[i] - B12 *(Y2s[i] - A2ma[J1s[i],]) , A12[J1s[i],], S12[J1s[i],]) # here we use A2ma not A2s!
        lnorm4 = lognormpdf(Y4s[i] - B43 *(Y3s[i] - A3ma[J1s[i],]) , A43[J1s[i],], S43[J1s[i],]) # here we use A3ma not A3s!
        lall = ltau + lnorm4 + lnorm1 + lnorm2 + lnorm3
        if (any(is.na(lall))) { cat("lall is na!\n"); stop =T; break }
        liks = liks + logsumexp(lall)
        taus[i,] = exp(lall - logsumexp(lall))
      }

      # compute prior
      lik_prior = (dprior-1) * sum(log(pk0))
      lik = liks + lik_prior

    } else {
      cat("skiping first max step, using supplied posterior probabilities\n")
    }

    tic("estep")

    if (stop) break;

    # ---------- MAX STEP -------------#
    rws     = c(t(taus)) + ctrl$posterior_reg
    wv2     = c(t(S2s[J1s,]))^2
    wv32    = c(t(S3s[J1s,]))^2
    rwav    = c(rws/wv2,rws/wv32)

    update_levels = (step%%2)==0

    if (update_levels==TRUE) {

      XX     = rBind(
        cBind(  Dkj1,       Dkj0   ),
        cBind(  -B32s*Dkj1  , Dkj1 )
      )

      DY2s    = as.matrix(kronecker(Y2s             ,rep(1,nk)))
      DY3s    = as.matrix(kronecker(Y3s - B32s*Y2s  ,rep(1,nk)))
      DYY     = rBind(DY2s,DY3s)

      fit     = slm.wfitc(XX,DYY,rwav,CS)$solution
      A2s     = t(rdim(fit[1:(nk*nf)],nk,nf))
      A3s     = t(rdim(fit[(nk*nf+1):(2*nk*nf)],nk,nf))
      if (ctrl$check_lik) m4.mixt.check.stayers(Y1s,Y2s,Y3s,Y4s,J1s,nk,Ns,model1,A3s=A3s,A2s=A2s);

      # maximize the variances
      fitv    = slm.wfit(XXV,(DYY - XX%*%fit)^2,rwav)
      S2s     = sqrt(t(rdim(coefficients(fitv)[1:(nk*nf)],nk,nf)))
      S3s     = sqrt(t(rdim(coefficients(fitv)[(nk*nf+1):(2*nk*nf)],nk,nf)))
      if (ctrl$check_lik) m4.mixt.check.stayers(Y1s,Y2s,Y3s,Y4s,J1s,nk,Ns,model1,S2s=S2s,S3s=S3s);
      S2s[S2s<ctrl$sd_floor]=ctrl$sd_floor # having a variance of exacvtly 0 creates problem in the likelihood
      S3s[S3s<ctrl$sd_floor]=ctrl$sd_floor

    } else {

      # update the covariance
      fit     = lm.wfit( as.matrix(c(t(Y2s - A2s[J1s,]))), c(t(Y3s - A3s[J1s,])),rws/wv32)
      B32s    = coefficients(fit)
      if (ctrl$check_lik) m4.mixt.check.stayers(Y1s,Y2s,Y3s,Y4s,J1s,nk,Ns,model1,B32s=B32s);
    }

    tic("mstep-ols")

    # -------- PK probabilities ------------ #
    # --- stayers --- #
    for (l1 in 1:nf) {
      I = which(J1s == l1)
      if (length(I)>1) {
        pk0[l1,] = colSums(taus[I,])
      } else {
        pk0[l1,] = taus[I,]
      }
      pk0[l1,] = (pk0[l1,] + dprior-1)/sum(pk0[l1,]+dprior-1)
    }
    if (ctrl$check_lik) m4.mixt.check.stayers(Y1s,Y2s,Y3s,Y4s,J1s,nk,Ns,model1,pk0=pk0);

    tic("mstep-pks")

    # checking model fit
    if ((!any(is.na(model0))) & ((step %% ctrl$nplot) == (ctrl$nplot-1))) {
      rr = addmom(A2s, model0$A2s,  "A2s")
      rr = addmom(A3s, model0$A3s,  "A3s", rr)
      rr = addmom(S12, model0$S12,  "S12", rr,type="var")
      rr = addmom(S3s, model0$S3s,  "S3s", rr,type="var")
      rr = addmom(S2s, model0$S2s,  "S2s", rr,type="var")
      rr = addmom(S43, model0$S43,  "S43", rr,type="var")
      rr = addmom(pk0, model0$pk0,  "pk0", rr,type="pr")
      rr = addmom(c(model$B12,model$B32m,model$B43,B32s), c(model0$B12,model0$B32m,model0$B43,model0$B32s), "rhos", rr,type="rho")

      print(ggplot(rr,aes(x=val1,y=val2,color=type)) + geom_point() + facet_wrap(~name,scale="free") + theme_bw() + geom_abline(linetype=2))
    } else {
      if ((step %% ctrl$nplot) == (ctrl$nplot-1)) {
        mm = list(A12=A12,pk0=pk0,A2s=A2s,A3s=A3s)
        #plot.endo.stayers(mm)
      }
    }

    ## -------- check convergence ------- #
    dlik = (lik - lik_old)/abs(lik_old)
    lik_old = lik
    lik_best = pmax(lik_best,lik)
    if ( (step %% ctrl$ncat) == 0) flog.info("[%3i][%s] lik=%4.4f dlik=%4.4e (sm:%4.4e) likm=%4.4e rho=%4.4e",step,ctrl$textapp,lik,dlik,dlik_smooth,model$likm,B32s)
    if (step>10) {
      dlik_smooth = 0.5*dlik_smooth + 0.5*dlik
      if (abs(dlik_smooth)<ctrl$tol) break;
    }

    tic("loop-wrap")
  }
  flog.info("[%3i][%s][final] lik=%4.4f dlik=%4.4e (sm:%4.4e) likm=%4.4e rho=%4.4e",step,ctrl$textapp,lik,dlik,dlik_smooth,model$likm,B32s)

  ## --stayers --
  model$pk0  = pk0
  model$A2s  = A2s
  model$S2s  = S2s
  model$A3s  = A3s
  model$S3s  = S3s
  model$B32s = B32s
  model$liks = liks
  model$NNs   = sdatae[,.N,j1][order(j1)][,N]

  end.time <- Sys.time()
  time.taken <- end.time - start.time

  return(list(tic = tic(), model=model,lik=lik,step=step,dlik=dlik,time.taken=time.taken,ctrl=ctrl,liks=liks,likm=likm))
}

#' The goal of this function is to recover rho1/rho2 from the
#' stayers.
#' @export
m4.mixt.rhoext.stayers.unc <- function(jdatae,sdatae,model,ctrl) {

  start.time <- Sys.time()
  tic <- tic.new()

  dprior = ctrl$dprior
  model0 = ctrl$model0
  tau = ctrl$tau

  ## ----- GET MODEL ---- #
  nk   = model$nk
  nf   = model$nf
  # Y1 | Y2
  A12  = model$A12
  B12  = model$B12
  S12  = model$S12
  A43  = model$A43
  B43  = model$B43
  S43  = model$S43
  ## --stayers --
  pk0  = model$pk0
  A2s  = model$A2s
  S2s  = model$S2s
  A3s  = model$A3s
  S3s  = model$S3s
  B32s = model$B32s

  ## ----- GET DATA ---- #
  # stayers
  Y1s = sdatae$y1
  Y2s = sdatae$y2
  Y3s = sdatae$y3
  Y4s = sdatae$y4
  J1s = sdatae$j1
  Ns  = sdatae[,.N]

  ## ----- PREPARE ---- #
  ## --- stayers --- #
  Dj1 = array(0,c(Ns,nf))
  Dj1[1:Ns + Ns*(J1s-1)]=1
  # duplicate them nk times (to use posterior weigths)
  Dkj1    = kronecker(Dj1,diag(nk))
  Dkj1    = as.matrix.csr(Dkj1)
  Dkj0    = Dkj1*0

  # get the constraints
  CS12 = cons.pad(cons.get(ctrl$cstr_type[1],ctrl$cstr_val[1],nk,nf),nk*nf*0, nk*nf*3)
  CS2  = cons.pad(cons.get(ctrl$cstr_type[2],ctrl$cstr_val[2],nk,nf),nk*nf*1, nk*nf*2)
  CS32 = cons.pad(cons.get(ctrl$cstr_type[3],ctrl$cstr_val[3],nk,nf),nk*nf*2, nk*nf*1)
  CS43 = cons.pad(cons.get(ctrl$cstr_type[4],ctrl$cstr_val[4],nk,nf),nk*nf*3, nk*nf*0)
  # combine them
  CS = cons.bind(cons.bind(cons.bind(CS12,CS2),CS32),CS43)

  # create the stationary contraints
  if (ctrl$fixb==T) {
    CS2 = cons.fixb(nk,nf,4)
    CS  = cons.bind(CS2,CS)
  }

  # create a constraint for the variances
  if (ctrl$model_var==T) {
    CSw = cons.none(nk,nf*4)
  } else{
    CS12 = cons.pad(cons.mono_k(nk,nf),nk*nf*0, nk*nf*3)
    CS2  = cons.pad(cons.mono_k(nk,nf),nk*nf*1, nk*nf*2)
    CS32 = cons.pad(cons.mono_k(nk,nf),nk*nf*2, nk*nf*1)
    CS43 = cons.pad(cons.mono_k(nk,nf),nk*nf*3, nk*nf*0)
    CSw  = cons.bind(cons.bind(cons.bind(CS12,CS2),CS32),CS43)
    CSw$meq = length(CSw$H)
  }

  # ---- prepare regressions covariates ---- #

  XXV     = rBind(
      cBind(  Dkj1, Dkj0, Dkj0, Dkj0 ),
      cBind(  Dkj0, Dkj1, Dkj0, Dkj0 ),
      cBind(  Dkj0, Dkj0, Dkj1, Dkj0 ),
      cBind(  Dkj0, Dkj0, Dkj0, Dkj1 )
  )

  lik_old  = -Inf
  lik      = -Inf
  lik_best = -Inf
  likm=0
  dlik_ma = 1

  model$dprior = dprior
  model1 = copy(model)

  tic("prep")

  stop = F;
  for (step in 1:ctrl$maxiter) {


    model1[c('A2s','A3s','S3s','S2s','pk0','B32s','A12','A43','S12','S43','B12','B43')] = list(A2s,A3s,S3s,S2s,pk0,B32s,A12,A43,S12,S43,B12,B43)

    # ---------- E STEP ------------- #
    # compute the tau probabilities and the likelihood
    if (is.na(tau[1]) | (step>1)) {

      # --- stayers --- #
      taus = array(0,c(Ns,nk))
      liks = 0
      for (i in 1:Ns) {
        ltau = log(pk0[J1s[i],])
        lnorm2 = lognormpdf(Y2s[i]                                , A2s[J1s[i],], S2s[J1s[i],])
        lnorm3 = lognormpdf(Y3s[i] - B32s*(Y2s[i] - A2s[J1s[i],]) , A3s[J1s[i],], S3s[J1s[i],])
        lnorm1 = lognormpdf(Y1s[i] - B12 *(Y2s[i] - A2s[J1s[i],]) , A12[J1s[i],], S12[J1s[i],])
        lnorm4 = lognormpdf(Y4s[i] - B43 *(Y3s[i] - A3s[J1s[i],]) , A43[J1s[i],], S43[J1s[i],])
        lall = ltau + lnorm4 + lnorm1 + lnorm2 + lnorm3
        if (any(is.na(lall))) { cat("lall is na!\n"); stop =T; break }
        liks = liks + logsumexp(lall)
        taus[i,] = exp(lall - logsumexp(lall))
      }

      # compute prior
      lik_prior = (dprior-1) * sum(log(pk0))
      lik = liks + lik_prior

    } else {
      cat("skiping first max step, using supplied posterior probabilities\n")
    }

    tic("estep")

    if (stop) break;

    # ---------- MAX STEP ------------- #
    rws     = c(t(taus+0.00000))
    wv12    = c(t(S12[J1s,]))^2
    wv2     = c(t(S2s[J1s,]))^2
    wv32    = c(t(S3s[J1s,]))^2
    wv43    = c(t(S43[J1s,]))^2
    rwav    = c(rws/wv12,rws/wv2,rws/wv32,rws/wv43)

    update_levels = ((step%%2)==0) | (step<5)

    if (update_levels==TRUE) {

     XX     = rBind(
        cBind(  Dkj1, -B12*Dkj1  , Dkj0      , Dkj0 ),
        cBind(  Dkj0,      Dkj1  , Dkj0      , Dkj0 ),
        cBind(  Dkj0, -B32s*Dkj1 , Dkj1      , Dkj0 ),
        cBind(  Dkj0,      Dkj0  , -B43*Dkj1 , Dkj1 )
      )

      DY1s    = as.matrix(kronecker(Y1s - B12 *Y2s  ,rep(1,nk)))
      DY2s    = as.matrix(kronecker(Y2s             ,rep(1,nk)))
      DY3s    = as.matrix(kronecker(Y3s - B32s*Y2s  ,rep(1,nk)))
      DY4s    = as.matrix(kronecker(Y4s - B43 *Y3s  ,rep(1,nk)))
      DYY     = rBind(DY1s,DY2s,DY3s,DY4s)

      fit     = slm.wfitc(XX,DYY,rwav,CS)$solution
      A12     = t(rdim(fit[1:(nk*nf)],nk,nf))
      A2s     = t(rdim(fit[(nk*nf+1):(2*nk*nf)],nk,nf))
      A3s     = t(rdim(fit[(2*nk*nf+1):(3*nk*nf)],nk,nf))
      A43     = t(rdim(fit[(3*nk*nf+1):(4*nk*nf)],nk,nf))
      if (ctrl$check_lik) m4.mixt.check.stayers.unc(Y1s,Y2s,Y3s,Y4s,J1s,nk,Ns,model1,A3s=A3s,A2s=A2s,A12=A12,A43=A43);
      lres = rdim(DYY - XX%*%fit,nk,Ns,4)

      # maximize the variances
      fitv  = slm.wfitc(XXV,(DYY - XX%*%fit)^2,rwav,CSw)$solution
      S12   = sqrt(t(rdim(fitv[    1      :(  nk*nf)],nk,nf)))
      S2s   = sqrt(t(rdim(fitv[(  nk*nf+1):(2*nk*nf)],nk,nf)))
      S3s   = sqrt(t(rdim(fitv[(2*nk*nf+1):(3*nk*nf)],nk,nf)))
      S43   = sqrt(t(rdim(fitv[(3*nk*nf+1):(4*nk*nf)],nk,nf)))
      if (ctrl$check_lik) m4.mixt.check.stayers.unc(Y1s,Y2s,Y3s,Y4s,J1s,nk,Ns,model1,S12=S12,S3s=S3s,S2s=S2s,S43=S43);

    } else {

      # update the covariance
      fit     = lm.wfit( as.matrix(c(t(Y2s - A2s[J1s,]))), c(t(Y3s - A3s[J1s,])),rws/wv32)
      B32s    = coefficients(fit)
      if (ctrl$check_lik) m4.mixt.check.stayers.unc(Y1s,Y2s,Y3s,Y4s,J1s,nk,Ns,model1,B32s=B32s);

      if (ctrl$est_rho[2]==TRUE) {
        fit     = lm.wfit( as.matrix(c(t(Y3s - A3s[J1s,]))), c(t(Y4s - A43[J1s,])),rws/wv43)
        B43     = coefficients(fit)
        if (ctrl$check_lik) m4.mixt.check.stayers.unc(Y1s,Y2s,Y3s,Y4s,J1s,nk,Ns,model1,B43=B43);
      }

      if (ctrl$est_rho[1]==TRUE) {
        fit     = lm.wfit( as.matrix(c(t(Y2s - A2s[J1s,]))), c(t(Y1s - A12[J1s,])),rws/wv12)
        B12     = coefficients(fit)
        if (ctrl$check_lik) m4.mixt.check.stayers.unc(Y1s,Y2s,Y3s,Y4s,J1s,nk,Ns,model1,B12=B12);
      }
    }

    tic("mstep-ols")

    # -------- PK probabilities ------------ #
    # --- stayers ---
    for (l1 in 1:nf) {
      I = which(J1s == l1)
      if (length(I)>1) {
        pk0[l1,] = colSums(taus[I,])
      } else {
        pk0[l1,] = taus[I,]
      }
      pk0[l1,] = (pk0[l1,] + dprior-1)/sum(pk0[l1,]+dprior-1)
    }
    if (ctrl$check_lik) m4.mixt.check.stayers.unc(Y1s,Y2s,Y3s,Y4s,J1s,nk,Ns,model1,pk0=pk0);

    tic("mstep-pks")

    # checking model fit
    if ((!any(is.na(model0))) & ((step %% ctrl$nplot) == (ctrl$nplot-1))) {
      rr = addmom(A2s, model0$A2s,  "A2s")
      rr = addmom(A3s, model0$A3s,  "A3s", rr)
      rr = addmom(S12, model0$S12,  "S12", rr,type="var")
      rr = addmom(S3s, model0$S3s,  "S3s", rr,type="var")
      rr = addmom(S2s, model0$S2s,  "S2s", rr,type="var")
      rr = addmom(S43, model0$S43,  "S43", rr,type="var")
      rr = addmom(pk0, model0$pk0,  "pk0", rr,type="pr")
      rr = addmom(c(model$B12,model$B32m,model$B43,B32s), c(model0$B12,model0$B32m,model0$B43,model0$B32s), "rhos", rr,type="rho")

      print(ggplot(rr,aes(x=val1,y=val2,color=type)) + geom_point() + facet_wrap(~name,scale="free") + theme_bw() + geom_abline(linetype=2))
    } else {
      if ((step %% ctrl$nplot) == (ctrl$nplot-1)) {
        mm = list(A12=A12,pk0=pk0,A2s=A2s,A3s=A3s)
        plot.endo.stayers(mm)
      }
    }

    # -------- check convergence ------- #
    dlik = (lik - lik_old)/abs(lik_old)
    lik_old = lik
    lik_best = pmax(lik_best,lik)
    dlik_smooth = 0
    catf("[%3i] lik=%4.4f dlik=%4.4e (%4.4e) dlik0=%4.4e liks=%4.4e likm=%4.4e r12=%.2f r32=%.2f r43=%.2f\n",step,lik,dlik,dlik_smooth,lik_best-lik,liks,model$likm,B12,B32s,B43)
    if (step>10) {
      dlik_smooth = 0.5*dlik_smooth + 0.5*dlik
      if (abs(dlik_smooth)<ctrl$tol) break;
    }

    tic("loop-wrap")
  }

  # compute some moments
  #rrd = data.frame()
  rrd2 = data.frame()
  ww  = rdim(rws,nk,Ns);
  for (iie in 1:4) {
    for (iik in 1:nk) {
      for (iil in 1:nf) {
        XS = scale(sample(lres[iik,J1s==iil,iie],sum(J1s==iil),prob=ww[iik,J1s==iil],replace=TRUE))
        #qs = quantile(XS,prob=seq(0.001,0.999,l=100))
        #qn = qnorm(seq(0.001,0.999,l=100))
        #rrd = rbind(rrd,data.frame(qs=qs,qn=qn,k=iik,l=iil,eq=iie))
        rrd2 = rbind(rrd2,data.frame(m=mean(XS),v=var(XS),s=skewness(XS),kurt=kurtosis(XS), k=iik,l=iil,eq=iie))
  }}}
  #  browser()
  #  ggplot(subset(rrd,eq==2),aes(x=qn,y=qs)) + geom_line()+geom_point() + geom_abline(linetype=2) + facet_grid(k~l,scales="free")+theme_bw()



  ## --stayers --
  model$liks = liks
  model[c('A2s','A3s','S3s','S2s','pk0','B32s','A12','A43','S12','S43','B12','B43')] = list(A2s,A3s,S3s,S2s,pk0,B32s,A12,A43,S12,S43,B12,B43)

  end.time <- Sys.time()
  time.taken <- end.time - start.time

  return(list(tic = tic(), model=model,lik=lik,step=step,tau=tau,dlik=dlik,time.taken=time.taken,ctrl=ctrl,liks=liks,likm=likm,res_distr=rrd2))
}



#' The goal of this function is to recover rho1/rho2 from the
#' stayers. This estimates non normal error distributions.
#' @export
m4.mixt.rhoext.stayers.unc.np <- function(jdatae,sdatae,model,ctrl) {

  start.time <- Sys.time()
  tic <- tic.new()

  dprior = ctrl$dprior
  model0 = ctrl$model0
  tau = ctrl$tau

  # ----- GET MODEL ---- #
  nk   = model$nk
  nf   = model$nf
  # Y1 | Y2
  A12  = model$A12
  B12  = model$B12
  S12  = model$S12
  A43  = model$A43
  B43  = model$B43
  S43  = model$S43
  ## --stayers --
  pk0  = model$pk0
  A2s  = model$A2s
  S2s  = model$S2s
  A3s  = model$A3s
  S3s  = model$S3s
  B32s = model$B32s

  # ----- GET DATA ---- #
  # stayers
  Y1s = sdatae$y1
  Y2s = sdatae$y2
  Y3s = sdatae$y3
  Y4s = sdatae$y4
  J1s = sdatae$j1
  Ns  = sdatae[,.N]

  # ----- PREPARE ---- #
  # --- stayers --- #
  Dj1 = array(0,c(Ns,nf))
  Dj1[1:Ns + Ns*(J1s-1)]=1
  # duplicate them nk times (to use posterior weigths)
  Dkj1    = kronecker(Dj1,diag(nk))
  Dkj1    = as.matrix.csr(Dkj1)
  Dkj0    = Dkj1*0

  # get the constraints
  CS12 = cons.pad(cons.get(ctrl$cstr_type[1],ctrl$cstr_val[1],nk,nf),nk*nf*0, nk*nf*3)
  CS2  = cons.pad(cons.get(ctrl$cstr_type[2],ctrl$cstr_val[2],nk,nf),nk*nf*1, nk*nf*2)
  CS32 = cons.pad(cons.get(ctrl$cstr_type[3],ctrl$cstr_val[3],nk,nf),nk*nf*2, nk*nf*1)
  CS43 = cons.pad(cons.get(ctrl$cstr_type[4],ctrl$cstr_val[4],nk,nf),nk*nf*3, nk*nf*0)
  # combine them
  CS = cons.bind(cons.bind(cons.bind(CS12,CS2),CS32),CS43)

  # create the stationary contraints
  if (ctrl$fixb==T) {
    CS2 = cons.fixb(nk,nf,4)
    CS  = cons.bind(CS2,CS)
  }

  # create a constraint for the variances
  if (ctrl$model_var==T) {
    CSw = cons.none(nk,nf*4)
  } else{
    CS12 = cons.pad(cons.mono_k(nk,nf),nk*nf*0, nk*nf*3)
    CS2  = cons.pad(cons.mono_k(nk,nf),nk*nf*1, nk*nf*2)
    CS32 = cons.pad(cons.mono_k(nk,nf),nk*nf*2, nk*nf*1)
    CS43 = cons.pad(cons.mono_k(nk,nf),nk*nf*3, nk*nf*0)
    CSw  = cons.bind(cons.bind(cons.bind(CS12,CS2),CS32),CS43)
    CSw$meq = length(CSw$H)
  }

  # ---- prepare regressions covariates ---- #

  XXV     = rBind(
      cBind(  Dkj1, Dkj0, Dkj0, Dkj0 ),
      cBind(  Dkj0, Dkj1, Dkj0, Dkj0 ),
      cBind(  Dkj0, Dkj0, Dkj1, Dkj0 ),
      cBind(  Dkj0, Dkj0, Dkj0, Dkj1 )
  )

  # ---- prepare the the error densities ---- #
  # we start with normal(0,1)
  nkd = 64
  XR = rnorm(10000)
  kd = density(XR,n=nkd)
  apprs = list()
  for (t in 1:4) {
    apprs[[t]] = approxfun(kd$x,kd$y,rule=2)
  }
  logkernpdf <- function(z,mu,sigma,t) log(  apprs[[t]]((z-mu)/sigma) )
  kernpdf <- function(z,mu,sigma,t) (  apprs[[t]]((z-mu)/sigma) )

  lik_old  = -Inf
  lik      = -Inf
  lik_best = -Inf
  likm=0
  dlik_ma = 1

  model$dprior = dprior
  model1 = copy(model)

  tic("prep")

  stop = F;
  for (step in 1:ctrl$maxiter) {


    model1[c('A2s','A3s','S3s','S2s','pk0','B32s','A12','A43','S12','S43','B12','B43')] = list(A2s,A3s,S3s,S2s,pk0,B32s,A12,A43,S12,S43,B12,B43)

    # ---------- E STEP ------------- #
    # compute the tau probabilities and the likelihood
    if (is.na(tau[1]) | (step>1)) {

      # --- stayers --- #
      taus = array(0,c(Ns,nk))
      liks = 0
      for (i in 1:Ns) {
        ltau = log(pk0[J1s[i],])
        # lnorm2 = lognormpdf(Y2s[i]                                , A2s[J1s[i],], S2s[J1s[i],])
        # lnorm3 = lognormpdf(Y3s[i] - B32s*(Y2s[i] - A2s[J1s[i],]) , A3s[J1s[i],], S3s[J1s[i],])
        # lnorm1 = lognormpdf(Y1s[i] - B12 *(Y2s[i] - A2s[J1s[i],]) , A12[J1s[i],], S12[J1s[i],])
        # lnorm4 = lognormpdf(Y4s[i] - B43 *(Y3s[i] - A3s[J1s[i],]) , A43[J1s[i],], S43[J1s[i],])
        lnorm2 = logkernpdf(Y2s[i]                                , A2s[J1s[i],], S2s[J1s[i],],2)
        lnorm3 = logkernpdf(Y3s[i] - B32s*(Y2s[i] - A2s[J1s[i],]) , A3s[J1s[i],], S3s[J1s[i],],3)
        lnorm1 = logkernpdf(Y1s[i] - B12 *(Y2s[i] - A2s[J1s[i],]) , A12[J1s[i],], S12[J1s[i],],1)
        lnorm4 = logkernpdf(Y4s[i] - B43 *(Y3s[i] - A3s[J1s[i],]) , A43[J1s[i],], S43[J1s[i],],4)
        lall = ltau + lnorm4 + lnorm1 + lnorm2 + lnorm3
        if (any(!is.finite(lall))) { cat("lall is na!\n"); browser(); stop =T; break }
        liks = liks + logsumexp(lall)
        taus[i,] = exp(lall - logsumexp(lall))
      }

      # compute prior
      lik_prior = (dprior-1) * sum(log(pk0))
      lik = liks + lik_prior

    } else {
      cat("skiping first max step, using supplied posterior probabilities\n")
    }

    tic("estep")

    if (stop) break;

    # ---------- MAX STEP ------------- #
    rws     = c(t(taus+0.00000))
    wv12    = c(t(S12[J1s,]))^2
    wv2     = c(t(S2s[J1s,]))^2
    wv32    = c(t(S3s[J1s,]))^2
    wv43    = c(t(S43[J1s,]))^2
    rwav    = c(rws/wv12,rws/wv2,rws/wv32,rws/wv43)

    update_levels = ((step%%2)==0) | (step<5)

    if (update_levels==TRUE) {

     XX     = rBind(
        cBind(  Dkj1, -B12*Dkj1  , Dkj0      , Dkj0 ),
        cBind(  Dkj0,      Dkj1  , Dkj0      , Dkj0 ),
        cBind(  Dkj0, -B32s*Dkj1 , Dkj1      , Dkj0 ),
        cBind(  Dkj0,      Dkj0  , -B43*Dkj1 , Dkj1 )
      )

      DY1s    = as.matrix(kronecker(Y1s - B12 *Y2s  ,rep(1,nk)))
      DY2s    = as.matrix(kronecker(Y2s             ,rep(1,nk)))
      DY3s    = as.matrix(kronecker(Y3s - B32s*Y2s  ,rep(1,nk)))
      DY4s    = as.matrix(kronecker(Y4s - B43 *Y3s  ,rep(1,nk)))
      DYY     = rBind(DY1s,DY2s,DY3s,DY4s)

      #browser()
      fit     = slm.wfitc(XX,DYY,rwav,CS)$solution
      A12     = t(rdim(fit[1:(nk*nf)],nk,nf))
      A2s     = t(rdim(fit[(nk*nf+1):(2*nk*nf)],nk,nf))
      A3s     = t(rdim(fit[(2*nk*nf+1):(3*nk*nf)],nk,nf))
      A43     = t(rdim(fit[(3*nk*nf+1):(4*nk*nf)],nk,nf))
      if (ctrl$check_lik) m4.mixt.check.stayers.unc(Y1s,Y2s,Y3s,Y4s,J1s,nk,Ns,model1,A3s=A3s,A2s=A2s,A12=A12,A43=A43);

      # maximize the variances
      fitv  = slm.wfitc(XXV,(DYY - XX%*%fit)^2,rwav,CSw)$solution
      S12   = sqrt(t(rdim(fitv[    1      :(  nk*nf)],nk,nf)))
      S2s   = sqrt(t(rdim(fitv[(  nk*nf+1):(2*nk*nf)],nk,nf)))
      S3s   = sqrt(t(rdim(fitv[(2*nk*nf+1):(3*nk*nf)],nk,nf)))
      S43   = sqrt(t(rdim(fitv[(3*nk*nf+1):(4*nk*nf)],nk,nf)))
      if (ctrl$check_lik) m4.mixt.check.stayers.unc(Y1s,Y2s,Y3s,Y4s,J1s,nk,Ns,model1,S12=S12,S3s=S3s,S2s=S2s,S43=S43);

      # update the marginal distributions ( we construct residuals for each time period with correct weights)
      browser()
      RS = DYY - XX%*%fit
      Nt = length(rws)
      for (it in 1:4) {
        kd1 = density(RS[((it-1)*Nt + 1):(it*Nt)],weights = rws/sum(rws))
        apprs[[it]] = approxfun(kd1$x,kd1$y,rule=2)
      }

    } else {

      # update the covariance
      fit     = lm.wfit( as.matrix(c(t(Y2s - A2s[J1s,]))), c(t(Y3s - A3s[J1s,])),rws/wv32)
      B32s    = coefficients(fit)
      if (ctrl$check_lik) m4.mixt.check.stayers.unc(Y1s,Y2s,Y3s,Y4s,J1s,nk,Ns,model1,B32s=B32s);

      if (ctrl$est_rho[2]==TRUE) {
        fit     = lm.wfit( as.matrix(c(t(Y3s - A3s[J1s,]))), c(t(Y4s - A43[J1s,])),rws/wv43)
        B43     = coefficients(fit)
        if (ctrl$check_lik) m4.mixt.check.stayers.unc(Y1s,Y2s,Y3s,Y4s,J1s,nk,Ns,model1,B43=B43);
      }

      if (ctrl$est_rho[1]==TRUE) {
        fit     = lm.wfit( as.matrix(c(t(Y2s - A2s[J1s,]))), c(t(Y1s - A12[J1s,])),rws/wv12)
        B12     = coefficients(fit)
        if (ctrl$check_lik) m4.mixt.check.stayers.unc(Y1s,Y2s,Y3s,Y4s,J1s,nk,Ns,model1,B12=B12);
      }
    }

    tic("mstep-ols")

    # -------- PK probabilities ------------ #
    # --- stayers --- #
    for (l1 in 1:nf) {
      I = which(J1s == l1)
      if (length(I)>1) {
        pk0[l1,] = colSums(taus[I,])
      } else {
        pk0[l1,] = taus[I,]
      }
      pk0[l1,] = (pk0[l1,] + dprior-1)/sum(pk0[l1,]+dprior-1)
    }
    if (ctrl$check_lik) m4.mixt.check.stayers.unc(Y1s,Y2s,Y3s,Y4s,J1s,nk,Ns,model1,pk0=pk0);

    tic("mstep-pks")

    # checking model fit
    if ((!any(is.na(model0))) & ((step %% ctrl$nplot) == (ctrl$nplot-1))) {
      rr = addmom(A2s, model0$A2s,  "A2s")
      rr = addmom(A3s, model0$A3s,  "A3s", rr)
      rr = addmom(S12, model0$S12,  "S12", rr,type="var")
      rr = addmom(S3s, model0$S3s,  "S3s", rr,type="var")
      rr = addmom(S2s, model0$S2s,  "S2s", rr,type="var")
      rr = addmom(S43, model0$S43,  "S43", rr,type="var")
      rr = addmom(pk0, model0$pk0,  "pk0", rr,type="pr")
      rr = addmom(c(model$B12,model$B32m,model$B43,B32s), c(model0$B12,model0$B32m,model0$B43,model0$B32s), "rhos", rr,type="rho")

      print(ggplot(rr,aes(x=val1,y=val2,color=type)) + geom_point() + facet_wrap(~name,scale="free") + theme_bw() + geom_abline(linetype=2))
    } else {
      if ((step %% ctrl$nplot) == (ctrl$nplot-1)) {
        mm = list(A12=A12,pk0=pk0,A2s=A2s,A3s=A3s)
        plot.endo.stayers(mm)
      }
    }

    # -------- check convergence ------- #
    dlik = (lik - lik_old)/abs(lik_old)
    lik_old = lik
    lik_best = pmax(lik_best,lik)
    dlik_smooth = 0
    catf("[%3i] lik=%4.4f dlik=%4.4e (%4.4e) dlik0=%4.4e liks=%4.4e likm=%4.4e r12=%.2f r32=%.2f r43=%.2f\n",step,lik,dlik,dlik_smooth,lik_best-lik,liks,model$likm,B12,B32s,B43)
    if (step>10) {
      dlik_smooth = 0.5*dlik_smooth + 0.5*dlik
      if (abs(dlik_smooth)<ctrl$tol) break;
    }

    tic("loop-wrap")
  }

  ## --stayers --
  model$liks = liks
  model[c('A2s','A3s','S3s','S2s','pk0','B32s','A12','A43','S12','S43','B12','B43')] = list(A2s,A3s,S3s,S2s,pk0,B32s,A12,A43,S12,S43,B12,B43)

  end.time <- Sys.time()
  time.taken <- end.time - start.time

  return(list(tic = tic(), model=model,lik=lik,step=step,tau=tau,dlik=dlik,time.taken=time.taken,ctrl=ctrl,liks=liks,likm=likm))
}





#' @export
em.end.level.shock <- function(model,l=0) {

  shockMat <- function(A,l) {

    if (typeof(A)=="double") {
      S = runif(length(A))
      A = A * (1 + 2*l*(S-.5))
    } else {
      tsize = cumprod(dim(A))
      tsize = tsize[length(tsize)]
      S = array(runif(tsize),dim(A))
      A = A * (1 + 2*l*(S-.5))
    }
    return(A)
  }

  model = within(model, {
    A12  = shockMat(A12,l)
    S12  = shockMat(S12,l)
    A43  = shockMat(A43,l)
    S43  = shockMat(S43,l)
    A2ma = shockMat(A2ma,l)
    A3ma = shockMat(A3ma,l)
    S2m  = shockMat(S2m,l)
    S3m  = shockMat(S3m,l)

    A3mb = shockMat(A3mb,l)
    A2mb = shockMat(A2mb,l)

    B32m = shockMat(B32m,l)
    B43  = shockMat(B43,l)
    B12  = shockMat(B12,l)

    pk1  = shockMat(pk1,l)
    pk0  = shockMat(pk1,l)
  })

  return(model)
}


#' @export
m4.mixt.movers.plotw <-function(model) {
  g1 = wplot(model$A12) + ggtitle("A12")
  g2 = wplot(model$A2ma) + ggtitle("A2ma")
  g3 = wplot(model$A3ma) + ggtitle("A3ma")
  g4 = wplot(model$A43) + ggtitle("A43")
  multiplot(g1,g2,g3,g4,cols = 2)
}

#' @export
m4.mixt.stayers.plotw <-function(model) {
  g1 = wplot(model$A2s) + ggtitle("A2s")
  g2 = wplot(model$A3s) + ggtitle("A3s")
  g3 = wplot(model$A12) + ggtitle("A12")
  g4 = pplot(model$pk0) + ggtitle("pk0")
  multiplot(g1,g2,g3,g4,cols = 2)
}

#' @export
endo.compare <- function(model1,model2) {

  rr = addmom(model1$A2s, model2$A2s,    "A2s")
  rr = addmom(model1$A3s, model2$A3s,    "A3s", rr)
  rr = addmom(model1$A12, model2$A12,    "A12", rr)
  rr = addmom(model1$A43, model2$A43,    "A43", rr)
  rr = addmom(model1$A2ma, model2$A2ma,  "A2ma", rr)
  rr = addmom(model1$A2mb, model2$A2mb,  "A2mb", rr)
  rr = addmom(model1$A3ma, model2$A3ma,  "A3ma", rr)
  rr = addmom(model1$A3mb, model2$A3mb,  "A3mb", rr)
  rr = addmom(model1$pk0, model2$pk0,  "pk0", rr,type="pr")
  rr = addmom(model1$pk1, model2$pk1,  "pk1", rr,type="pr")
  rr = addmom(c(model1$B12,model1$B32m,model1$B43,model1$B32s), c(model2$B12,model2$B32m,model2$B43,model2$B32s), "rhos", rr,type="rho")

  print(ggplot(rr,aes(x=val1,y=val2,color=type)) + geom_point() +
          facet_wrap(~name,scale="free") + theme_bw() + geom_abline(linetype=2)+ xlab("model1") + ylab("model2"))

  catf("model1 likm= %4.4e \t liks= %4.4e \n",model1$likm,model1$liks)
  catf("model2 likm= %4.4e \t liks= %4.4e \n",model2$likm,model2$liks)

}

test.neglik <- function() {

}

# testing function getMean
test.getMean <- function() {

  # --- movers ---
  computeLikSmall <- function(Y2m,A2ma,A2mb,S2m,pk1,dprior,getTau=F) {
    taum = array(0,c(Nm,nk))
    likm = 0
    for (i in 1:Nm) {
      ltau = log(pk1[JJm[i],])
      lnorm2 = lognormpdf(Y2m[i]              , A2ma[J1m[i],] + A2mb[J2m[i]], S2m[J1m[i],])
      lnorm3 = lognormpdf(Y3m[i] - B32m*Y2m[i], A3ma[J2m[i],] + A3mb[J1m[i]], S3m[J2m[i],])
      lnorm1 = lognormpdf(Y1m[i] - B12 *Y2m[i], A12[J1m[i],], S12[J1m[i],])
      lnorm4 = lognormpdf(Y4m[i] - B43 *Y3m[i], A43[J2m[i],], S43[J2m[i],])
      lall = ltau + lnorm2 +lnorm3
      likm = likm + logsumexp(lall)
      taum[i,] = exp(lall - logsumexp(lall))
    }

    # compute prior
    lik_prior = (dprior-1) * sum(log(pk1))
    lik = likm + lik_prior
    if (getTau) return(taum);
    return(lik)
  }

  rwm = c(t(computeLikSmall(Y2m,A2ma,A2mb,S2m,pk1,dprior,getTau=T)))
  # getMean finds the minimum least square solution with weights
  # under some constraints. The padding should jsut treat the last variable
  # without constraints.
  fit = getMean(Dkj1_j2,DY2m,rwm,C1,H1,meq,nk,nf,pad=nf-1)

  computeLikSmall(Y2m,A2ma,A2mb,S2m,pk1,dprior)
  computeLikSmall(Y2m,fit$A,A2mb,S2m,pk1,dprior)
  computeLikSmall(Y2m,A2ma,c(0,fit$B),S2m,pk1,dprior)


}

#' Estimates the dynamic model parameters, but keeps the rhos fixed
#'
#' @export
m4.mixtr.rhoext.movers <- function(jdatae,sdatae,model,ctrl,em.order=1) {

  start.time <- Sys.time()
  tic <- tic.new()

  dprior = ctrl$dprior
  model0 = ctrl$model0
  taum = ctrl$tau

  ### ----- GET MODEL  ---
  nk   = model$nk
  nf   = model$nf
  # Y1 | Y2
  A12  = model$A12
  B12  = model$B12
  S12  = model$S12
  A43  = model$A43
  B43  = model$B43
  S43  = model$S43
  ## --movers
  pk1  = model$pk1
  A2ma = model$A2ma
  A2mb = model$A2mb
  S2m  = model$S2m
  A3ma = model$A3ma
  A3mb = model$A3mb
  A2ma = model$A2ma
  A2mb = model$A2mb
  S3m  = model$S3m
  B32m = model$B32m

  # ----- GET DATA
  # movers
  Y1m = jdatae$y1
  Y2m = jdatae$y2
  Y3m = jdatae$y3
  Y4m = jdatae$y4
  J1m = jdatae$j1
  J2m = jdatae$j2
  JJm = J1m + nf*(J2m-1)
  Nm = jdatae[,.N]

  # Transform pk1 and use NNm
  dtmp = jdatae[,.N,list(j1,j2)]
  dtmp = dtmp[,list(pr=pk1[j1 + nf*(j2-1),],k=1:nk,jj=j1 + nf*(j2-1),N),list(j1,j2)]
  dtmp[,V1:=pr*N/sum(pr*N),k]
  pk1  = acast(dtmp, jj ~ k,value.var = "V1")
  pkk = dtmp[,sum(pr*N),k][order(k),V1]
  pkk = pkk/sum(pkk)

  # get the constraints
  CS12 = cons.pad(cons.get(ctrl$cstr_type[1],ctrl$cstr_val[1],nk,nf),nk*nf*0, nk*nf*3 + (nf-1)*2)
  CS2  = cons.pad(cons.get(ctrl$cstr_type[2],ctrl$cstr_val[2],nk,nf),nk*nf*1, nk*nf*2 + (nf-1)*2)
  CS32 = cons.pad(cons.get(ctrl$cstr_type[3],ctrl$cstr_val[3],nk,nf),nk*nf*2, nk*nf*1 + (nf-1)*2)
  CS43 = cons.pad(cons.get(ctrl$cstr_type[4],ctrl$cstr_val[4],nk,nf),nk*nf*3, nk*nf*0 + (nf-1)*2)
  # combine them
  CS = cons.bind(cons.bind(cons.bind(CS12,CS2),CS32),CS43)

  # create the stationary contraints
  if (ctrl$fixb==T) {
    CS2 = cons.pad(cons.fixb(nk,nf,4),0,(nf-1)*2)
    CS  = cons.bind(CS2,CS)
  }

  # create a constraint for the variances
  if (ctrl$model_var==T) {
    CSw = cons.none(nk,nf*4)
  } else{
    CS12 = cons.pad(cons.mono_k(nk,nf),nk*nf*0, nk*nf*3)
    CS2  = cons.pad(cons.mono_k(nk,nf),nk*nf*1, nk*nf*2)
    CS32 = cons.pad(cons.mono_k(nk,nf),nk*nf*2, nk*nf*1)
    CS43 = cons.pad(cons.mono_k(nk,nf),nk*nf*3, nk*nf*0)
    CSw  = cons.bind(cons.bind(cons.bind(CS12,CS2),CS32),CS43)
    CSw$meq = length(CSw$H)
  }

  # prepare matrices aggregated at the type level
  Dkj1f     = diag(nf)  %x% rep(1,nf) %x% diag(nk)            # A[k,l] coefficients for j1
  Dkj2f     = rep(1,nf) %x% diag(nf)  %x% diag(nk)            # A[k,l] coefficients for j2
  Dj1f_bis  = (diag(nf)  %x% rep(1,nf) %x% rep(1,nk))[,2:nf]  # C[l]   coefficients for j1
  Dj2f_bis  = (rep(1,nf) %x% diag(nf)  %x% rep(1,nk))[,2:nf]  # C[l]   coefficients for j2

  # regression matrix for the variance
  XXV = rBind(
    cBind(    Dkj1f, 0*Dkj1f,  0*Dkj2f, 0*Dkj2f ),
    cBind(  0*Dkj1f,   Dkj1f,  0*Dkj2f, 0*Dkj2f ),
    cBind(  0*Dkj1f, 0*Dkj1f,    Dkj2f, 0*Dkj2f ),
    cBind(  0*Dkj1f, 0*Dkj1f,  0*Dkj2f,   Dkj2f )
  )

  ## --- prepare regressions covariates --- #
  # create the depend variables

  lik_old  = -Inf
  lik      = -Inf
  lik_best = -Inf
  liks = 0
  likm=0

  lpt1 = array(0,c(Nm,nk))
  lpt2 = array(0,c(Nm,nk))
  lpt3 = array(0,c(Nm,nk))
  lpt4 = array(0,c(Nm,nk))
  lp   = array(0,c(Nm,nk))


  tic("prep")
  stop = F;
  for (step in 1:ctrl$maxiter) {

    model1 = list(nk=nk,nf=nk,A12=A12,B12=B12,S12=S12,
                  A2ma=A2ma,A2mb=A2mb,S2m=S2m,
                  A3ma=A3ma,A3mb=A3mb,S3m=S3m,B32m=B32m,
                  A43=A43,B43=B43,S43=S43,
                  pk1=pk1,dprior=dprior)

    ### ---------- E STEP ------------- #
    # compute the tau probabilities and the likelihood
    if (is.na(taum[1]) | (step>1)) {

      # for efficiency we want to group by (l1,l2)
      for (l1 in 1:nf) for (l2 in 1:nf) {
        I = which( (J1m==l1) & (J2m==l2))
        ll = l1 + nf*(l2-1)
        if (length(I)==0) next;

        for (k in 1:nk) {
          lpt2[I,k] = lognormpdf(Y2m[I]                                        , A2ma[l1,k] + A2mb[l2] , S2m[l1,k])
          lpt1[I,k] = lognormpdf(Y1m[I] - B12 *(Y2m[I] - A2ma[l1,k])           , A12[l1,k], S12[l1,k])
          lpt3[I,k] = lognormpdf(Y3m[I] - B32m*(Y2m[I] - A2ma[l1,k] - A2mb[l2]), A3ma[l2,k] + A3mb[l1], S3m[l2,k])
          lpt4[I,k] = lognormpdf(Y4m[I] - B43 *(Y3m[I] - A3ma[l2,k])           , A43[l2,k], S43[l2,k])

          # sum the log of the periods
          lp[I,k] = log(pkk[k]) + log(pk1[ll,k]) + lpt1[I,k] + lpt2[I,k] + lpt3[I,k] + lpt4[I,k]
        }
      }

      liks     = sum(logRowSumExp(lp))
      taum     = exp(lp - spread(logRowSumExp(lp),2,nk)) # normalize the k probabilities Pr(k|Y1,Y2,Y3,Y4,l)

      # compute prior for both p(j1,j2|l) and p(l)
      lik_prior = (dprior[1]-1) * sum(log(pk1)) + (dprior[2]-1)*sum(log(pkk))
      lik = liks + lik_prior

    } else {
      cat("skiping first max step, using supplied posterior probabilities\n")
    }

    tic("estep")

    if (stop) break;

    # ---------- MAX STEP ------------- #
    # taum = makePosteriorStochastic(tau = taum,m = ctrl$stochastic) # if we want to implement stochastic EM

    # we start by recovering the posterior weight, and the variances for each term
    rwm  = c(t(taum + ctrl$posterior_reg))
    wv12 = c(t(S12[J1m,]))^2
    wv2  = c(t(S2m[J1m,]))^2
    wv32 = c(t(S3m[J2m,]))^2
    wv43 = c(t(S43[J2m,]))^2

    # here we switch between updating the rho paramerters
    # and updating the other "level" parameters
    update_levels = ((em.order + (step))%%2 )==0
    if (all(ctrl$est_rho[1:3]==FALSE)) update_levels=TRUE;

    if (update_levels==TRUE) {
      DYY  = array(0,c(nk,nf,nf,4))
      WWT  = array(1e-7,c(nk,nf,nf,4))

      for (l1 in 1:nf) for (l2 in 1:nf) {
        I = which( (J1m==l1) & (J2m==l2))
        if (length(I)==0) next;
        for (k in 1:nk) {
          # compute the posterior weight, it's not time specific
          ww = sum(taum[I,k])

          # construct dependent for each time period k,l2,l1,
          DYY[k,l2,l1,1] = sum(  (Y1m[I] - B12  * Y2m[I]) * taum[I,k] )/ww
          DYY[k,l2,l1,2] = sum(  (Y2m[I]                ) * taum[I,k] )/ww
          DYY[k,l2,l1,3] = sum(  (Y3m[I] - B32m * Y2m[I]) * taum[I,k] )/ww
          DYY[k,l2,l1,4] = sum(  (Y4m[I] - B43  * Y3m[I]) * taum[I,k] )/ww

          # Scaling the weight by the time specific variance
          WWT[k,l2,l1,1] = ww/S12[l1,k]^2
          WWT[k,l2,l1,2] = ww/S2m[l1,k]^2
          WWT[k,l2,l1,3] = ww/S3m[l2,k]^2
          WWT[k,l2,l1,4] = ww/S43[l2,k]^2
        }
      }

      # modifiy the regressor to take into account the possibly changed B12, B23m, B43
      XXf     = rBind(
        cBind(    Dkj1f,  -B12*Dkj1f,       0*Dkj1f,      0*Dkj1f,      0*Dj2f_bis,  0*Dj1f_bis),
        cBind(  0*Dkj1f,       Dkj1f,       0*Dkj1f,      0*Dkj1f,        Dj2f_bis,  0*Dj1f_bis),
        cBind(  0*Dkj1f, -B32m*Dkj1f,         Dkj2f,      0*Dkj1f,  -B32m*Dj2f_bis,    Dj1f_bis),
        cBind(  0*Dkj1f,     0*Dkj1f,    -B43*Dkj2f,        Dkj2f,      0*Dj2f_bis,  0*Dj1f_bis)
      )

      if (any(!is.finite(WWT))) browser();
      if (any(WWT==0)) browser();

      #fit   = lm.wfitc(XXf,as.numeric(DYY),as.numeric(WWT),CS$C,CS$H,CS$meq)$solution
      fit   = slm.wfitc(XXf,as.numeric(DYY),as.numeric(WWT),CS)$solution
      A12   = t(rdim(fit[1:(nk*nf)],nk,nf))
      A2ma  = t(rdim(fit[(nk*nf+1):(2*nk*nf)],nk,nf))
      A3ma  = t(rdim(fit[(2*nk*nf+1):(3*nk*nf)],nk,nf))
      A43   = t(rdim(fit[(3*nk*nf+1):(4*nk*nf)],nk,nf))
      A2mb  = c(0,fit[(4*nk*nf+1):(4*nk*nf+nf-1)])
      A3mb  = c(0,fit[(4*nk*nf+nf):(4*nk*nf+2*(nf-1))])
      if (ctrl$check_lik) {
        m4.mixtr.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,A12=A12,A2ma=A2ma,A3ma=A3ma,A43=A43,A2mb=A2mb,A3mb=A3mb)
      }

      # compute the variances!!!!
      DYY_bar = array(0,c(nk,nf,nf,4))
      DYY_bar[] = XXf%*%fit
      DYYV   =  array(0,c(nk,nf,nf,4))

      for (l1 in 1:nf) for (l2 in 1:nf) {
        I = which( (J1m==l1) & (J2m==l2))
        if (length(I)==0) next;
        for (k in 1:nk) {
          # construct dependent for each time period k,l2,l1,
          ww = sum(taum[I,k])
          DYYV[k,l2,l1,1] = sum(  (Y1m[I] - B12  * Y2m[I] - DYY_bar[k,l2,l1,1])^2 * taum[I,k] )/ww
          DYYV[k,l2,l1,2] = sum(  (Y2m[I]                 - DYY_bar[k,l2,l1,2])^2 * taum[I,k] )/ww
          DYYV[k,l2,l1,3] = sum(  (Y3m[I] - B32m * Y2m[I] - DYY_bar[k,l2,l1,3])^2 * taum[I,k] )/ww
          DYYV[k,l2,l1,4] = sum(  (Y4m[I] - B43  * Y3m[I] - DYY_bar[k,l2,l1,4])^2 * taum[I,k] )/ww
        }
      }

      fitv  = slm.wfitc(XXV,as.numeric(DYYV),as.numeric(WWT),CSw)$solution
      S12   = sqrt(t(rdim(fitv[1:(nk*nf)],nk,nf)))
      S2m   = sqrt(t(rdim(fitv[(nk*nf+1):(2*nk*nf)],nk,nf)))
      S3m   = sqrt(t(rdim(fitv[(2*nk*nf+1):(3*nk*nf)],nk,nf)))
      S43   = sqrt(t(rdim(fitv[(3*nk*nf+1):(4*nk*nf)],nk,nf)))
      if (ctrl$check_lik) {
        m4.mixtr.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,A12=A12,A2ma=A2ma,A3ma=A3ma,A43=A43,A2mb=A2mb,A3mb=A3mb,S12=S12,S2m=S2m,S3m=S3m,S43=S43)
      }

    } else {
      # estimate the rhos -- use the other parameters
      if (ctrl$est_rho[1]==TRUE) {
        fit = lm.wfit( as.matrix(c(t(Y2m - A2ma[J1m,]))), c(t(Y1m - A12[J1m,])),rwm/wv12)
        B12 = coefficients(fit)
        if (ctrl$check_lik) { m4.mixtr.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,B12=B12)}
      }
      if (ctrl$est_rho[3]==TRUE) {
        fit = lm.wfit( as.matrix(c(t(Y3m - A3ma[J2m,]))), c(t(Y4m - A43[J2m,])),rwm/wv43)
        B43 = coefficients(fit)
        if (ctrl$check_lik) { m4.mixtr.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,B43=B43)}
      }
      if (ctrl$est_rho[2]==TRUE) {
        fit = lm.wfit( as.matrix(c(t(Y2m - A2ma[J1m,] - A2mb[J2m]))), c(t(Y3m - A3ma[J2m,] - A3mb[J1m])),rwm/wv32)
        B32m = coefficients(fit)
        if (ctrl$check_lik) { m4.mixtr.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,B32m=B32m)}
      }
    }

    tic("mstep-ols")

    ## -------- PK probabilities ------------ #
    ## --- movers --- #
    for (l1 in 1:nf) for (l2 in 1:nf) {
      jj = l1 + nf*(l2-1)
      I = which(JJm == jj)
      if (length(I)>1) {
        pk1[jj,] = colSums(taum[I,])
      } else if (length(I)==0) { # this deals with the case where the cell is empty
        pk1[jj,] = 1/nk
      } else {
        pk1[jj,] = taum[I,]
      }
    }
    for (k in 1:nk) {
      pkk[k]  = sum(pk1[,k])
      pk1[,k] = (pk1[,k] + dprior[1]-1 )/(sum(pk1[,k] + dprior[1] -1 ))
    }
    #pkk = pkk/sum(pkk)
    pkk = (pkk + dprior[2]-1 )/(sum(pkk + dprior[2] -1 ))
    if (ctrl$check_lik) { m4.mixtr.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,pk1=pk1)}

    tic("mstep-pks")

    # checking model fit
    if ((!any(is.na(model0))) & ((step %% ctrl$nplot) == (ctrl$nplot-1))) {
      I1 = order(colSums(A2ma))
      I2 = order(colSums(model0$A2ma))
      rr = addmom(A2ma[,I1],model0$A2ma[,I2],"A2ma")
      rr = addmom(A2mb,model0$A2mb,"A2mb",rr)
      rr = addmom(A3ma[,I1],model0$A3ma[,I2],"A3ma",rr)
      rr = addmom(A3mb,model0$A3mb,"A3mb",rr)
      rr = addmom(A12[,I1],model0$A12[,I2],"A12",rr)
      rr = addmom(A43[,I1],model0$A43[,I2],"A43",rr)
      rr = addmom(c(B12,B32m,B43), c(model0$B12,model0$B32m,model0$B43), "rhos", rr,type="rho")
      rr = addmom(S2m[,I1], model0$S2m[,I2], "S2m", rr,type="var")
      rr = addmom(S3m[,I1], model0$S3m[,I2], "S3m", rr,type="var")
      rr = addmom(S12[,I1],model0$S12[,I2],"S12",rr,type="var")
      rr = addmom(S43[,I1],model0$S43[,I2],"S43",rr,type="var")

      rr = addmom(pk1,model0$pk1,"pk1",rr,type="pr")

      print(ggplot(rr,aes(x=val2,y=val1,color=type)) + geom_point() + facet_wrap(~name,scale="free") + theme_bw() + geom_abline(linetype=2))
    } else {
      if ((step %% ctrl$nplot) == (ctrl$nplot-1)) {
        mm = list(A12=A12,A2ma=A2ma,A3ma=A3ma,A43=A43)
        plot.endo.movers(mm)
      }
    }

    # -------- check convergence ------- #
    dlik = (lik - lik_old)/abs(lik_old)
    lik_old = lik
    lik_best = pmax(lik_best,lik)
    if ( (step %% ctrl$ncat) == 0) {
      flog.info("[%3i] levels=%i lik=%4.4f dlik=%4.4e dlik0=%4.4e liks=%4.4e likm=%4.4e r12=%.3f r23=%.3f r43=%.3f",step,update_levels+0,lik,dlik,lik_best-lik,liks,likm,B12,B32m,B43);
      flog.info( paste(round(pkk,3),collapse = " "))
    }
    if (step>10) if (abs(dlik)<ctrl$tol) break;

    tic("loop-wrap")
  }

  # udpdate the pk1 to match usual definition
  model$pkbis_j1j2 = pk1
  model$pkbis_k    = pkk
  for (l1 in 1:nf)  for (l2 in 1:nf){
    ll = l1 + nf*(l2-1)
    pk1[ll,] = pk1[ll,] * pkk/sum(pk1[ll,] * pkk)
  }

  # Y1 | Y2
  model$A12  = A12
  model$B12  = B12
  model$S12  = S12
  model$A43  = A43
  model$B43  = B43
  model$S43  = S43
  ## --movers --
  model$pk1  = pk1
  model$A2ma = A2ma
  model$A2mb = A2mb
  model$S2m  = S2m
  model$A3ma = A3ma
  model$A3mb = A3mb
  model$S3m  = S3m
  model$B32m = B32m

  model$NNm = acast(jdatae[,.N,list(j1,j2)],j1~j2,fill=0,value.var="N")
  model$likm=lik

  end.time <- Sys.time()
  time.taken <- end.time - start.time

  return(list(tic = tic(), model=model,lik=lik,step=step,tau=taum,dlik=dlik,time.taken=time.taken,ctrl=ctrl,liks=liks,likm=likm))
}


#' Estimates the dynamic model parameters, but keeps the rhos fixed.
#' This uses a parametrization p(k'|k,l) and p(k|l)
#'
#' @export
m4.mixtr2.rhoext.movers <- function(jdatae,sdatae,model,ctrl,em.order=1) {

  start.time <- Sys.time()
  tic <- tic.new()

  dprior = ctrl$dprior
  model0 = ctrl$model0
  taum = ctrl$tau

  ### ----- GET MODEL  ---
  nk   = model$nk
  nf   = model$nf
  # Y1 | Y2
  A12  = model$A12
  B12  = model$B12
  S12  = model$S12
  A43  = model$A43
  B43  = model$B43
  S43  = model$S43
  ## --movers
  pk1  = model$pk1
  A2ma = model$A2ma
  A2mb = model$A2mb
  S2m  = model$S2m
  A3ma = model$A3ma
  A3mb = model$A3mb
  A2ma = model$A2ma
  A2mb = model$A2mb
  S3m  = model$S3m
  B32m = model$B32m

  # ----- GET DATA
  # movers
  Y1m = jdatae$y1
  Y2m = jdatae$y2
  Y3m = jdatae$y3
  Y4m = jdatae$y4
  J1m = jdatae$j1
  J2m = jdatae$j2
  JJm = J1m + nf*(J2m-1)
  Nm = jdatae[,.N]

  # Transform pk1 and use NNm
  dtmp = jdatae[,.N,list(j1,j2)]
  dtmp = dtmp[,list(pr=pk1[j1 + nf*(j2-1),],k=1:nk,jj=j1 + nf*(j2-1),N),list(j1,j2)]


  # get Pr[j1,k], only if it is not in the model
  if ("pk_j1k" %in% names(model)) {
    pk_j1k = model$pk_j1k
  } else {
    pk_j1k  = acast(dtmp[,sum(pr*N),list(j1,k)], j1 ~ k,value.var = "V1")
  }

  # get Pr[j2|k,j]
  if ("pk_j2_j1k" %in% names(model)) {
    pk_j2_j1k = model$pk_j2_j1k
  } else {
    dtmp[,V1:=pr*N/sum(pr*N),list(j1,k)]
    pk_j2_j1k  = acast(dtmp, j2 ~ j1 ~ k,value.var = "V1")
  }

  # construct prior for Pr[j2|k,j]
  dtmp = jdatae[,.N,list(j1,j2)]
  alpha_j1j2 = acast(dtmp,j1~j2,value.var = "N")

  # remove pk1 to make sure we don't use it anywhere
  rm(pk1)

  # get the constraints
  CS12 = cons.pad(cons.get(ctrl$cstr_type[1],ctrl$cstr_val[1],nk,nf),nk*nf*0, nk*nf*3 + (nf-1)*2)
  CS2  = cons.pad(cons.get(ctrl$cstr_type[2],ctrl$cstr_val[2],nk,nf),nk*nf*1, nk*nf*2 + (nf-1)*2)
  CS32 = cons.pad(cons.get(ctrl$cstr_type[3],ctrl$cstr_val[3],nk,nf),nk*nf*2, nk*nf*1 + (nf-1)*2)
  CS43 = cons.pad(cons.get(ctrl$cstr_type[4],ctrl$cstr_val[4],nk,nf),nk*nf*3, nk*nf*0 + (nf-1)*2)
  # combine them
  CS = cons.bind(cons.bind(cons.bind(CS12,CS2),CS32),CS43)

  # create the stationary contraints
  if (ctrl$fixb==T) {
    CS2 = cons.pad(cons.fixb(nk,nf,4),0,(nf-1)*2)
    CS  = cons.bind(CS2,CS)
  }

  # create a constraint for the variances
  if (ctrl$model_var==T) {
    CSw = cons.none(nk,nf*4)
  } else{
    CS12 = cons.pad(cons.mono_k(nk,nf),nk*nf*0, nk*nf*3)
    CS2  = cons.pad(cons.mono_k(nk,nf),nk*nf*1, nk*nf*2)
    CS32 = cons.pad(cons.mono_k(nk,nf),nk*nf*2, nk*nf*1)
    CS43 = cons.pad(cons.mono_k(nk,nf),nk*nf*3, nk*nf*0)
    CSw  = cons.bind(cons.bind(cons.bind(CS12,CS2),CS32),CS43)
    CSw$meq = length(CSw$H)
  }

  # prepare matrices aggregated at the type level
  Dkj1f     = diag(nf)  %x% rep(1,nf) %x% diag(nk)            # A[k,l] coefficients for j1
  Dkj2f     = rep(1,nf) %x% diag(nf)  %x% diag(nk)            # A[k,l] coefficients for j2
  Dj1f_bis  = (diag(nf)  %x% rep(1,nf) %x% rep(1,nk))[,2:nf]  # C[l]   coefficients for j1
  Dj2f_bis  = (rep(1,nf) %x% diag(nf)  %x% rep(1,nk))[,2:nf]  # C[l]   coefficients for j2

  # regression matrix for the variance
  XXV = rBind(
    cBind(    Dkj1f, 0*Dkj1f,  0*Dkj2f, 0*Dkj2f ),
    cBind(  0*Dkj1f,   Dkj1f,  0*Dkj2f, 0*Dkj2f ),
    cBind(  0*Dkj1f, 0*Dkj1f,    Dkj2f, 0*Dkj2f ),
    cBind(  0*Dkj1f, 0*Dkj1f,  0*Dkj2f,   Dkj2f )
  )

  ## --- prepare regressions covariates --- #
  # create the depend variables

  lik_old  = -Inf
  lik      = -Inf
  lik_best = -Inf
  liks = 0
  likm=0

  lpt1 = array(0,c(Nm,nk))
  lpt2 = array(0,c(Nm,nk))
  lpt3 = array(0,c(Nm,nk))
  lpt4 = array(0,c(Nm,nk))
  lp   = array(0,c(Nm,nk))


  tic("prep")
  stop = F;
  for (step in 1:ctrl$maxiter) {

    model1 = list(nk=nk,nf=nk,A12=A12,B12=B12,S12=S12,
                  A2ma=A2ma,A2mb=A2mb,S2m=S2m,
                  A3ma=A3ma,A3mb=A3mb,S3m=S3m,B32m=B32m,
                  A43=A43,B43=B43,S43=S43,
                  dprior=dprior)

    ### ---------- E STEP ------------- #
    # compute the tau probabilities and the likelihood
    if (is.na(taum[1]) | (step>1)) {

      lik_prior_pkk = 0;

      # for efficiency we want to group by (l1,l2)
      for (l1 in 1:nf) for (l2 in 1:nf) {
        I = which( (J1m==l1) & (J2m==l2))
        ll = l1 + nf*(l2-1)
        if (length(I)==0) next;

        for (k in 1:nk) {
          lpt2[I,k] = lognormpdf(Y2m[I]                                        , A2ma[l1,k] + A2mb[l2] , S2m[l1,k])
          lpt1[I,k] = lognormpdf(Y1m[I] - B12 *(Y2m[I] - A2ma[l1,k])           , A12[l1,k], S12[l1,k])
          lpt3[I,k] = lognormpdf(Y3m[I] - B32m*(Y2m[I] - A2ma[l1,k] - A2mb[l2]), A3ma[l2,k] + A3mb[l1], S3m[l2,k])
          lpt4[I,k] = lognormpdf(Y4m[I] - B43 *(Y3m[I] - A3ma[l2,k])           , A43[l2,k], S43[l2,k])

          # sum the log of the periods
          lp[I,k] = log(pk_j1k[l1,k]) + log(pk_j2_j1k[l2,l1,k]) + lpt1[I,k] + lpt2[I,k] + lpt3[I,k] + lpt4[I,k]
        }

        # adding the prior on Pr[j2|j1,k]
        lik_prior_pkk = lik_prior_pkk +
          (dprior[2]*alpha_j1j2[l1,l2]-1)*sum(log(pk_j2_j1k[l2,l1,]))
      }

      liks     = sum(logRowSumExp(lp))
      taum     = exp(lp - spread(logRowSumExp(lp),2,nk)) # normalize the k probabilities Pr(k|Y1,Y2,Y3,Y4,l)

      # adding the prior on Pr[j1,k]
      lik_prior_pkk = lik_prior_pkk + (dprior[1]-1) * sum(log(pk_j1k))
      lik = liks + lik_prior_pkk

    } else {
      cat("skiping first max step, using supplied posterior probabilities\n")
    }

    tic("estep")

    if (stop) break;

    # ---------- MAX STEP ------------- #
    # taum = makePosteriorStochastic(tau = taum,m = ctrl$stochastic) # if we want to implement stochastic EM

    # we start by recovering the posterior weight, and the variances for each term
    rwm  = c(t(taum + ctrl$posterior_reg))
    wv12 = c(t(S12[J1m,]))^2
    wv2  = c(t(S2m[J1m,]))^2
    wv32 = c(t(S3m[J2m,]))^2
    wv43 = c(t(S43[J2m,]))^2

    # here we switch between updating the rho paramerters
    # and updating the other "level" parameters
    update_levels = ((em.order + (step))%%2 )==0
    if (all(ctrl$est_rho[1:3]==FALSE)) update_levels=TRUE;

    if (update_levels==TRUE) {
      DYY  = array(0,c(nk,nf,nf,4))
      WWT  = array(1e-7,c(nk,nf,nf,4))

      for (l1 in 1:nf) for (l2 in 1:nf) {
        I = which( (J1m==l1) & (J2m==l2))
        if (length(I)==0) next;
        for (k in 1:nk) {
          # compute the posterior weight, it's not time specific
          ww = sum(taum[I,k])

          # construct dependent for each time period k,l2,l1,
          DYY[k,l2,l1,1] = sum(  (Y1m[I] - B12  * Y2m[I]) * taum[I,k] )/ww
          DYY[k,l2,l1,2] = sum(  (Y2m[I]                ) * taum[I,k] )/ww
          DYY[k,l2,l1,3] = sum(  (Y3m[I] - B32m * Y2m[I]) * taum[I,k] )/ww
          DYY[k,l2,l1,4] = sum(  (Y4m[I] - B43  * Y3m[I]) * taum[I,k] )/ww

          # Scaling the weight by the time specific variance
          WWT[k,l2,l1,1] = ww/S12[l1,k]^2
          WWT[k,l2,l1,2] = ww/S2m[l1,k]^2
          WWT[k,l2,l1,3] = ww/S3m[l2,k]^2
          WWT[k,l2,l1,4] = ww/S43[l2,k]^2
        }
      }

      # modifiy the regressor to take into account the possibly changed B12, B23m, B43
      XXf     = rBind(
        cBind(    Dkj1f,  -B12*Dkj1f,       0*Dkj1f,      0*Dkj1f,      0*Dj2f_bis,  0*Dj1f_bis),
        cBind(  0*Dkj1f,       Dkj1f,       0*Dkj1f,      0*Dkj1f,        Dj2f_bis,  0*Dj1f_bis),
        cBind(  0*Dkj1f, -B32m*Dkj1f,         Dkj2f,      0*Dkj1f,  -B32m*Dj2f_bis,    Dj1f_bis),
        cBind(  0*Dkj1f,     0*Dkj1f,    -B43*Dkj2f,        Dkj2f,      0*Dj2f_bis,  0*Dj1f_bis)
      )

      if (any(!is.finite(WWT))) browser();
      if (any(WWT==0)) browser();

      #fit   = lm.wfitc(XXf,as.numeric(DYY),as.numeric(WWT),CS$C,CS$H,CS$meq)$solution
      fit   = slm.wfitc(XXf,as.numeric(DYY),as.numeric(WWT),CS)$solution
      A12   = t(rdim(fit[1:(nk*nf)],nk,nf))
      A2ma  = t(rdim(fit[(nk*nf+1):(2*nk*nf)],nk,nf))
      A3ma  = t(rdim(fit[(2*nk*nf+1):(3*nk*nf)],nk,nf))
      A43   = t(rdim(fit[(3*nk*nf+1):(4*nk*nf)],nk,nf))
      A2mb  = c(0,fit[(4*nk*nf+1):(4*nk*nf+nf-1)])
      A3mb  = c(0,fit[(4*nk*nf+nf):(4*nk*nf+2*(nf-1))])
      if (ctrl$check_lik) {
        m4.mixtr.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,A12=A12,A2ma=A2ma,A3ma=A3ma,A43=A43,A2mb=A2mb,A3mb=A3mb)
      }

      # compute the variances!!!!
      DYY_bar = array(0,c(nk,nf,nf,4))
      DYY_bar[] = XXf%*%fit
      DYYV   =  array(0,c(nk,nf,nf,4))

      for (l1 in 1:nf) for (l2 in 1:nf) {
        I = which( (J1m==l1) & (J2m==l2))
        if (length(I)==0) next;
        for (k in 1:nk) {
          # construct dependent for each time period k,l2,l1,
          ww = sum(taum[I,k])
          DYYV[k,l2,l1,1] = sum(  (Y1m[I] - B12  * Y2m[I] - DYY_bar[k,l2,l1,1])^2 * taum[I,k] )/ww
          DYYV[k,l2,l1,2] = sum(  (Y2m[I]                 - DYY_bar[k,l2,l1,2])^2 * taum[I,k] )/ww
          DYYV[k,l2,l1,3] = sum(  (Y3m[I] - B32m * Y2m[I] - DYY_bar[k,l2,l1,3])^2 * taum[I,k] )/ww
          DYYV[k,l2,l1,4] = sum(  (Y4m[I] - B43  * Y3m[I] - DYY_bar[k,l2,l1,4])^2 * taum[I,k] )/ww
        }
      }

      fitv  = slm.wfitc(XXV,as.numeric(DYYV),as.numeric(WWT),CSw)$solution
      S12   = sqrt(t(rdim(fitv[1:(nk*nf)],nk,nf)))
      S2m   = sqrt(t(rdim(fitv[(nk*nf+1):(2*nk*nf)],nk,nf)))
      S3m   = sqrt(t(rdim(fitv[(2*nk*nf+1):(3*nk*nf)],nk,nf)))
      S43   = sqrt(t(rdim(fitv[(3*nk*nf+1):(4*nk*nf)],nk,nf)))
      if (ctrl$check_lik) {
        m4.mixtr.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,A12=A12,A2ma=A2ma,A3ma=A3ma,A43=A43,A2mb=A2mb,A3mb=A3mb,S12=S12,S2m=S2m,S3m=S3m,S43=S43)
      }

    } else {
      # estimate the rhos -- use the other parameters
      if (ctrl$est_rho[1]==TRUE) {
        fit = lm.wfit( as.matrix(c(t(Y2m - A2ma[J1m,]))), c(t(Y1m - A12[J1m,])),rwm/wv12)
        B12 = coefficients(fit)
        if (ctrl$check_lik) { m4.mixtr.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,B12=B12)}
      }
      if (ctrl$est_rho[3]==TRUE) {
        fit = lm.wfit( as.matrix(c(t(Y3m - A3ma[J2m,]))), c(t(Y4m - A43[J2m,])),rwm/wv43)
        B43 = coefficients(fit)
        if (ctrl$check_lik) { m4.mixtr.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,B43=B43)}
      }
      if (ctrl$est_rho[2]==TRUE) {
        fit = lm.wfit( as.matrix(c(t(Y2m - A2ma[J1m,] - A2mb[J2m]))), c(t(Y3m - A3ma[J2m,] - A3mb[J1m])),rwm/wv32)
        B32m = coefficients(fit)
        if (ctrl$check_lik) { m4.mixtr.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,B32m=B32m)}
      }
    }

    tic("mstep-ols")

    ## -------- PK probabilities ------------ #
    ## extract posterior probabilities Pr[j1,j2,k]
    pkk = array(0,c(nf,nf,nk))
    for (j1 in 1:nf) for (j2 in 1:nf) {
      jj = j1 + nf*(j2-1)
      I = which(JJm == jj)
      if (length(I)>1) {
        pkk[j1,j2,] = colSums(taum[I,])
      } else if (length(I)==0) { # this deals with the case where the cell is empty
        pkk[j1,j2,] = 1/nk
      } else {
        pkk[j1,j2,] = taum[I,]
      }
    }
    ## update P[j2|j1,k]
    for (j1 in 1:nf) for (k in 1:nk) {
      pk_j2_j1k[,j1,k]  = (pkk[j1,,k] + dprior[2]*alpha_j1j2[j1,]-1 )/(sum(pkk[j1,,k] + dprior[2]*alpha_j1j2[j1,]-1 ))
      pk_j1k[j1,k] = sum(pkk[j1,,k])
    }
    pk_j1k = (pk_j1k + dprior[1]-1 )/(sum( pk_j1k + dprior[1]-1  ))
    if (ctrl$check_lik) { m4.mixtr.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,pk1=pk1)}

    tic("mstep-pks")

    # checking model fit
    if ((!any(is.na(model0))) & ((step %% ctrl$nplot) == (ctrl$nplot-1))) {
      I1 = order(colSums(A2ma))
      I2 = order(colSums(model0$A2ma))
      rr = addmom(A2ma[,I1],model0$A2ma[,I2],"A2ma")
      rr = addmom(A2mb,model0$A2mb,"A2mb",rr)
      rr = addmom(A3ma[,I1],model0$A3ma[,I2],"A3ma",rr)
      rr = addmom(A3mb,model0$A3mb,"A3mb",rr)
      rr = addmom(A12[,I1],model0$A12[,I2],"A12",rr)
      rr = addmom(A43[,I1],model0$A43[,I2],"A43",rr)
      rr = addmom(c(B12,B32m,B43), c(model0$B12,model0$B32m,model0$B43), "rhos", rr,type="rho")
      rr = addmom(S2m[,I1], model0$S2m[,I2], "S2m", rr,type="var")
      rr = addmom(S3m[,I1], model0$S3m[,I2], "S3m", rr,type="var")
      rr = addmom(S12[,I1],model0$S12[,I2],"S12",rr,type="var")
      rr = addmom(S43[,I1],model0$S43[,I2],"S43",rr,type="var")

      rr = addmom(pk1,model0$pk1,"pk1",rr,type="pr")

      print(ggplot(rr,aes(x=val2,y=val1,color=type)) + geom_point() + facet_wrap(~name,scale="free") + theme_bw() + geom_abline(linetype=2))
    } else {
      if ((step %% ctrl$nplot) == (ctrl$nplot-1)) {
        mm = list(A12=A12,A2ma=A2ma,A3ma=A3ma,A43=A43)
        plot.endo.movers(mm)
      }
    }

    # -------- check convergence ------- #
    dlik = (lik - lik_old)/abs(lik_old)
    lik_old = lik
    lik_best = pmax(lik_best,lik)
    if ( (step %% ctrl$ncat) == 0) flog.info("[%3i] levels=%i lik=%4.4f dlik=%4.4e dlik0=%4.4e liks=%4.4e likm=%4.4e r12=%.3f r23=%.3f r43=%.3f",step,update_levels+0,lik,dlik,lik_best-lik,liks,likm,B12,B32m,B43);
    if (step>10) if (abs(dlik)<ctrl$tol) break;

    tic("loop-wrap")
  }

  # udpdate the pk1 to match usual definition
  model$pk_j2_j1k = pk_j2_j1k
  model$pk_j1k    = pk_j1k
  pk1 = array(0,c(nf*nf,nk))
  for (l1 in 1:nf)  for (l2 in 1:nf){
    ll = l1 + nf*(l2-1)
    pk1[ll,] = pk_j2_j1k[l2,l1,] * pk_j1k[l1,] / sum(pk_j2_j1k[l2,l1,] * pk_j1k[l1,])
  }

  # Y1 | Y2
  model$A12  = A12
  model$B12  = B12
  model$S12  = S12
  model$A43  = A43
  model$B43  = B43
  model$S43  = S43
  ## --movers --
  model$pk1  = pk1
  model$A2ma = A2ma
  model$A2mb = A2mb
  model$S2m  = S2m
  model$A3ma = A3ma
  model$A3mb = A3mb
  model$S3m  = S3m
  model$B32m = B32m

  model$NNm = acast(jdatae[,.N,list(j1,j2)],j1~j2,fill=0,value.var="N")
  model$likm=lik

  end.time <- Sys.time()
  time.taken <- end.time - start.time

  return(list(tic = tic(), model=model,lik=lik,step=step,tau=taum,dlik=dlik,time.taken=time.taken,ctrl=ctrl,liks=liks,likm=likm))
}


#' Estimates the dynamic model parameters, but keeps the rhos fixed. This uses a joint likelihood.
#'
#' @export
m4.mixt.rhoext.joint <- function(sim,model,ctrl,em.order=1) {

  jdatae  = sim$jdata
  sdatae  = sim$sdata

  start.time <- Sys.time()
  tic <- tic.new()

  dprior = ctrl$dprior
  model0 = ctrl$model0
  taum = ctrl$tau

  ### ----- GET MODEL  ---
  nk   = model$nk
  nf   = model$nf
  # Y1 | Y2
  A12  = model$A12
  B12  = model$B12
  S12  = model$S12
  A43  = model$A43
  B43  = model$B43
  S43  = model$S43
  ## --movers
  pk1  = model$pk1
  A2ma = model$A2ma
  A2mb = model$A2mb
  S2m  = model$S2m
  A3ma = model$A3ma
  A3mb = model$A3mb
  A2ma = model$A2ma
  A2mb = model$A2mb
  S3m  = model$S3m
  B32m = model$B32m
  ## -- stayers --
  pk0  = model$pk0
  A2s  = model$A2s
  S2s  = model$S2s
  A3s  = model$A3s
  S3s  = model$S3s
  B32s = model$B32s

  ## ----- GET DATA ---- #
  # stayers
  Y1s = sdatae$y1
  Y2s = sdatae$y2
  Y3s = sdatae$y3
  Y4s = sdatae$y4
  J1s = sdatae$j1
  Ns  = sdatae[,.N]
  # movers
  Y1m = jdatae$y1
  Y2m = jdatae$y2
  Y3m = jdatae$y3
  Y4m = jdatae$y4
  J1m = jdatae$j1
  J2m = jdatae$j2
  JJm = J1m + nf*(J2m-1)
  Nm = jdatae[,.N]


  # get the constraints for movers
  CS12  = cons.pad(cons.get(ctrl$cstr_type[1],ctrl$cstr_val[1],nk,nf),nk*nf*0, nk*nf*3 + (nf-1)*2 + 2*nk*nf)
  CS2   = cons.pad(cons.get(ctrl$cstr_type[2],ctrl$cstr_val[2],nk,nf),nk*nf*1, nk*nf*2 + (nf-1)*2 + 2*nk*nf)
  CS32  = cons.pad(cons.get(ctrl$cstr_type[3],ctrl$cstr_val[3],nk,nf),nk*nf*2, nk*nf*1 + (nf-1)*2 + 2*nk*nf)
  CS43  = cons.pad(cons.get(ctrl$cstr_type[4],ctrl$cstr_val[4],nk,nf),nk*nf*3, nk*nf*0 + (nf-1)*2 + 2*nk*nf)
  # get the constraints for stayers
  CS2s  = cons.pad(cons.get(ctrl$cstr_type[2],ctrl$cstr_val[2],nk,nf),nk*nf*4 + (nf-1)*2 , nk*nf)
  CS32s = cons.pad(cons.get(ctrl$cstr_type[3],ctrl$cstr_val[3],nk,nf),nk*nf*5 + (nf-1)*2 , 0)

  # combine them
  CS = CS12 %>%
      cons.bind(CS12) %>%
      cons.bind(CS2)  %>%
      cons.bind(CS32) %>%
      cons.bind(CS43) %>%
      cons.bind(CS2s) %>%
      cons.bind(CS32s)

  # create the stationary contraints
  if (ctrl$fixb==T) {
    CS2 = cons.pad(cons.fixb(nk,nf,4),0,(nf-1)*2+ 2*nk*nf)
    CS  = cons.bind(CS2,CS)
  }

  # create a constraint for the variances
  if (ctrl$model_var==T) {
    CSw = cons.none(nk,nf*6)
  } else{
    CS12 = cons.pad(cons.mono_k(nk,nf),nk*nf*0, nk*nf*5)
    CS2  = cons.pad(cons.mono_k(nk,nf),nk*nf*1, nk*nf*4)
    CS32 = cons.pad(cons.mono_k(nk,nf),nk*nf*2, nk*nf*3)
    CS43 = cons.pad(cons.mono_k(nk,nf),nk*nf*3, nk*nf*2)
    CS2s = cons.pad(cons.mono_k(nk,nf),nk*nf*4, nk*nf*1)
    CS3s = cons.pad(cons.mono_k(nk,nf),nk*nf*5, nk*nf*0)
    CSw  = cons.bind(cons.bind(cons.bind(cons.bind(cons.bind(CS12,CS2),CS32),CS43),CS2s),CS3s)
    CSw$meq = length(CSw$H)
  }

  # prepare regessors matrices aggregated at the type level for movers
  Dkj1f     = diag(nf)  %x% rep(1,nf) %x% diag(nk)            # A[k,l] coefficients for j1
  Dkj2f     = rep(1,nf) %x% diag(nf)  %x% diag(nk)            # A[k,l] coefficients for j2
  Dj1f_bis  = (diag(nf)  %x% rep(1,nf) %x% rep(1,nk))[,2:nf]  # C[l]   coefficients for j1
  Dj2f_bis  = (rep(1,nf) %x% diag(nf)  %x% rep(1,nk))[,2:nf]  # C[l]   coefficients for j2
  # prepare regessors matrices aggregated at the type level for stayers
  Dkj1s     = diag(nf)  %x% diag(nk)                          # A[k,l] coefficients for j1
  Dj1s_bis  = (diag(nf) %x% rep(1,nk))[,2:nf]                 # C[l] coefficients for j1

  # regression matrix for the variance
  XXV = rBind(
      cBind(    Dkj1f, 0*Dkj1f,  0*Dkj2f, 0*Dkj2f,  0*Dkj2f, 0*Dkj2f),
      cBind(  0*Dkj1f,   Dkj1f,  0*Dkj2f, 0*Dkj2f,  0*Dkj2f, 0*Dkj2f ),
      cBind(  0*Dkj1f, 0*Dkj1f,    Dkj2f, 0*Dkj2f,  0*Dkj2f, 0*Dkj2f ),
      cBind(  0*Dkj1f, 0*Dkj1f,  0*Dkj2f,   Dkj2f,  0*Dkj2f, 0*Dkj2f ),
      cBind(    Dkj1s, 0*Dkj1s,  0*Dkj1s, 0*Dkj1s,  0*Dkj1s, 0*Dkj1s ),  # I need to put this in diagonal! not stacked
      cBind(  0*Dkj1s, 0*Dkj1s,  0*Dkj1s, 0*Dkj1s,    Dkj1s, 0*Dkj1s ),
      cBind(  0*Dkj1s, 0*Dkj1s,  0*Dkj1s, 0*Dkj1s,  0*Dkj1s,   Dkj1s ),
      cBind(  0*Dkj1s, 0*Dkj1s,  0*Dkj1s,   Dkj1s,  0*Dkj1s, 0*Dkj1s )
  )

  ## --- prepare regressions covariates --- #
  # create the depend variables

  lik_old  = -Inf
  lik      = -Inf
  lik_best = -Inf
  liks=0
  likm=0

  # temporary variables to store posterior probabilities
  # for movers and stayers
  lpm   = array(0,c(Nm,nk))
  lps   = array(0,c(Ns,nk))

  tic("prep")
  stop = F;
  for (step in 1:ctrl$maxiter) {

    model1 = list(nk=nk,nf=nk,A12=A12,B12=B12,S12=S12,
                  A2ma=A2ma,A2mb=A2mb,S2m=S2m,
                  A3ma=A3ma,A3mb=A3mb,S3m=S3m,B32m=B32m,
                  A43=A43,B43=B43,S43=S43,
                  pk1=pk1,dprior=dprior)

    ### ---------- E STEP ------------- #

    # compute the tau probabilities and the likelihood for the movers
    if (is.na(taum[1]) | (step>1)) {

      #  -------- MOVERS ---------- #
      # for efficiency we want to group by (l1,l2)
      for (l1 in 1:nf) for (l2 in 1:nf) {
        I = which( (J1m==l1) & (J2m==l2))
        ll = l1 + nf*(l2-1)
        if (length(I)==0) next;

        for (k in 1:nk) {
          lpm_2 = lognormpdf(Y2m[I]                                        , A2ma[l1,k] + A2mb[l2] , S2m[l1,k])
          lpm_1 = lognormpdf(Y1m[I] - B12 *(Y2m[I] - A2ma[l1,k])           , A12[l1,k], S12[l1,k])
          lpm_3 = lognormpdf(Y3m[I] - B32m*(Y2m[I] - A2ma[l1,k] - A2mb[l2]), A3ma[l2,k] + A3mb[l1], S3m[l2,k])
          lpm_4 = lognormpdf(Y4m[I] - B43 *(Y3m[I] - A3ma[l2,k])           , A43[l2,k], S43[l2,k])

          # sum the log of the periods
          lpm[I,k] = log(pk1[ll,k]) + lpm_1 + lpm_2 + lpm_3 + lpm_4
        }
      }

      likm     = sum(logRowSumExp(lpm))
      taum     = exp(lpm - spread(logRowSumExp(lpm),2,nk)) # normalize the k probabilities Pr(k|Y1,Y2,Y3,Y4,l)

      # compute prior
      likm_prior = (dprior-1) * sum(log(pk1))
      likm = likm + likm_prior

      #  -------- STAYERS ---------- #
      taus = array(0,c(Ns,nk))
      liks = 0
      for (l1 in 1:nf) {
      	I = which(J1s==l1)

        for (k in 1:nk) {
          lps_2 = lognormpdf(Y2s[I]                                 , A2s[l1,k], S2s[l1,k])
          lps_1 = lognormpdf(Y1s[I] - B12 *(Y2s[I] -  A2s[l1,k])    , A12[l1,k], S12[l1,k])
          lps_3 = lognormpdf(Y3s[I] - B32s*(Y2s[I] -  A2s[l1,k])    , A3s[l1,k], S3s[l1,k])
          lps_4 = lognormpdf(Y4s[I] - B43 *(Y3s[I] -  A3s[l1,k])    , A43[l1,k], S43[l1,k])

          # sum the log of the periods
          lps[I,k] = log(pk0[l1,k]) + lps_1 + lps_2 + lps_3 + lps_4
        }
      }

      liks     = sum(logRowSumExp(lps))
      taus     = exp(lps - spread(logRowSumExp(lps),2,nk)) # normalize the k probabilities Pr(k|Y1,Y2,Y3,Y4,l)

      # compute prior
      liks_prior = (dprior-1) * sum(log(pk0))
      liks = liks + liks_prior
    } else {
      cat("skiping first max step, using supplied posterior probabilities\n")
    }

    lik_tot = ctrl$stayer_weight*liks + likm

    tic("estep")

    if (stop) break;

    # ---------- MAX STEP ------------- #
    # taum = makePosteriorStochastic(tau = taum,m = ctrl$stochastic) # if we want to implement stochastic EM


    # here we switch between updating the rho paramerters
    # and updating the other "level" parameters
    update_levels = ((em.order + (step))%%2 )==0
    if (all(ctrl$est_rho[1:3]==FALSE)) update_levels=TRUE;

    if (update_levels==TRUE) {
      DYYm  = array(0,c(nk,nf,nf,4))
      WWTm  = array(1e-7,c(nk,nf,nf,4))

      for (l1 in 1:nf) for (l2 in 1:nf) {
        I = which( (J1m==l1) & (J2m==l2))
        if (length(I)==0) next;
        for (k in 1:nk) {
          # compute the posterior weight, it's not time specific
          ww = sum(taum[I,k]) + 1e-7

          # construct dependent for each time period k,l2,l1,
          DYYm[k,l2,l1,1] = sum(  (Y1m[I] - B12  * Y2m[I]) * taum[I,k] )/ww
          DYYm[k,l2,l1,2] = sum(  (Y2m[I]                ) * taum[I,k] )/ww
          DYYm[k,l2,l1,3] = sum(  (Y3m[I] - B32m * Y2m[I]) * taum[I,k] )/ww
          DYYm[k,l2,l1,4] = sum(  (Y4m[I] - B43  * Y3m[I]) * taum[I,k] )/ww

          # Scaling the weight by the time specific variance
          WWTm[k,l2,l1,1] = ww/pmax(1e-7,S12[l1,k]^2)
          WWTm[k,l2,l1,2] = ww/pmax(1e-7,S2m[l1,k]^2)
          WWTm[k,l2,l1,3] = ww/pmax(1e-7,S3m[l2,k]^2)
          WWTm[k,l2,l1,4] = ww/pmax(1e-7,S43[l2,k]^2)
        }
      }

      DYYs  = array(0,c(nk,nf,4))
      WWTs  = array(1e-7,c(nk,nf,4))

      #for (k in 1:nk) taus[,k]=(sdatae$k==k); # @debug

      for (l1 in 1:nf) {
        I = which((J1s==l1))
        if (length(I)==0) next;
        for (k in 1:nk) {
          # compute the posterior weight, it's not time specific
          ww = sum(taus[I,k]) + 1e-7

          # construct dependent for each time period k,l2,l1,
          DYYs[k,l1,1] = sum(  (Y1s[I] - B12  * Y2s[I]) * taus[I,k] )/ww
          DYYs[k,l1,2] = sum(  (Y2s[I]                ) * taus[I,k] )/ww
          DYYs[k,l1,3] = sum(  (Y3s[I] - B32s * Y2s[I]) * taus[I,k] )/ww
          DYYs[k,l1,4] = sum(  (Y4s[I] - B43  * Y3s[I]) * taus[I,k] )/ww

          # Scaling the weight by the time specific variance
          WWTs[k,l1,1] = ww/pmax(1e-7,S12[l1,k]^2)
          WWTs[k,l1,2] = ww/pmax(1e-7,S2s[l1,k]^2)
          WWTs[k,l1,3] = ww/pmax(1e-7,S3s[l1,k]^2)
          WWTs[k,l1,4] = ww/pmax(1e-7,S43[l1,k]^2)
        }
      }

     # Constructing the full matrix of regressors.
     # the first 4 rows are for movers for each periods, the next 4 are for the stayers
     XXf     = rBind( #     A12,        A2ma,          A3ma,          A43,            A2mb,        A3mb,     A2s,     A3s
                cBind(    Dkj1f,  -B12*Dkj1f,       0*Dkj1f,      0*Dkj1f,      0*Dj2f_bis,  0*Dj1f_bis,     0*Dkj1f,     0*Dkj1f), # movers  T=1
                cBind(  0*Dkj1f,       Dkj1f,       0*Dkj1f,      0*Dkj1f,        Dj2f_bis,  0*Dj1f_bis,     0*Dkj1f,     0*Dkj1f), # movers  T=2
                cBind(  0*Dkj1f, -B32m*Dkj1f,         Dkj2f,      0*Dkj1f,  -B32m*Dj2f_bis,    Dj1f_bis,     0*Dkj1f,     0*Dkj1f), # movers  T=3
                cBind(  0*Dkj1f,     0*Dkj1f,    -B43*Dkj2f,        Dkj2f,      0*Dj2f_bis,  0*Dj1f_bis,     0*Dkj1f,     0*Dkj1f), # movers  T=4
                cBind(    Dkj1s,     0*Dkj1s,       0*Dkj1s,      0*Dkj1s,      0*Dj1s_bis,  0*Dj1s_bis,  -B12*Dkj1s,     0*Dkj1s), # movers  T=1
                cBind(  0*Dkj1s,     0*Dkj1s,       0*Dkj1s,      0*Dkj1s,      0*Dj1s_bis,  0*Dj1s_bis,       Dkj1s,     0*Dkj1s), # movers  T=2
                cBind(  0*Dkj1s,     0*Dkj1s,       0*Dkj1s,      0*Dkj1s,      0*Dj1s_bis,  0*Dj1s_bis, -B32s*Dkj1s,       Dkj1s), # movers  T=3
                cBind(  0*Dkj1s,     0*Dkj1s,       0*Dkj1s,        Dkj1s,      0*Dj1s_bis,  0*Dj1s_bis,     0*Dkj1s,  -B43*Dkj1s)  # movers  T=4
     )

     WW = c( as.numeric(WWTm), ctrl$stayer_weight * as.numeric(WWTs))
     fit   = slm.wfitc(XXf,c(as.numeric(DYYm),as.numeric(DYYs)),WW,CS)$solution
     is = 1
     A12   = t(rdim(fit[is:(is + nk*nf-1)],nk,nf)); is = is+nk*nf
     A2ma  = t(rdim(fit[is:(is + nk*nf-1)],nk,nf)); is = is+nk*nf
     A3ma  = t(rdim(fit[is:(is + nk*nf-1)],nk,nf)); is = is+nk*nf
     A43   = t(rdim(fit[is:(is + nk*nf-1)],nk,nf)); is = is+nk*nf
     A2mb  = c(0,fit[is:(is + nf-2)]);           is = is+nf-1
     A3mb  = c(0,fit[is:(is + nf-2)]);           is = is+nf-1
     A2s   = t(rdim(fit[is:(is + nk*nf-1)],nk,nf)); is = is+nk*nf
     A3s   = t(rdim(fit[is:(is + nk*nf-1)],nk,nf)); is = is+nk*nf
     if (ctrl$check_lik) {
       m4.mixt.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,A12=A12,A2ma=A2ma,A3ma=A3ma,A43=A43,A2mb=A2mb,A3mb=A3mb)
     }

    # variances - movers
    DYYm_bar   = array(0,c(nk,nf,nf,4))
    DYYm_bar[] = (XXf%*%fit)[1:(nk*nf*nf*4)]
    DYYVm      =  array(0,c(nk,nf,nf,4))
    for (l1 in 1:nf) for (l2 in 1:nf) {
      I = which( (J1m==l1) & (J2m==l2))
      if (length(I)==0) next;
      for (k in 1:nk) {
        # construct dependent for each time period k,l2,l1,
        ww = sum(taum[I,k]) + 1e-7
        DYYVm[k,l2,l1,1] = sum(  (Y1m[I] - B12  * Y2m[I] - DYYm_bar[k,l2,l1,1])^2 * taum[I,k] )/ww
        DYYVm[k,l2,l1,2] = sum(  (Y2m[I]                 - DYYm_bar[k,l2,l1,2])^2 * taum[I,k] )/ww
        DYYVm[k,l2,l1,3] = sum(  (Y3m[I] - B32m * Y2m[I] - DYYm_bar[k,l2,l1,3])^2 * taum[I,k] )/ww
        DYYVm[k,l2,l1,4] = sum(  (Y4m[I] - B43  * Y3m[I] - DYYm_bar[k,l2,l1,4])^2 * taum[I,k] )/ww
      }
    }

    # variances - stayers
    DYYs_bar   = array(0,c(nk,nf,4))
    DYYs_bar[] = (XXf%*%fit)[(nk*nf*nf*4+1):(nk*nf*nf*4 + nk*nf*4)]
    DYYVs      =  array(0,c(nk,nf,4))
    for (l1 in 1:nf) {
      I = which((J1s==l1))
      if (length(I)==0) next;
      for (k in 1:nk) {
        # construct dependent for each time period k,l2,l1,
        ww = sum(taus[I,k]) + 1e-7
        DYYVs[k,l1,1] = sum(  (Y1s[I] - B12  * Y2s[I] - DYYs_bar[k,l1,1])^2 * taus[I,k] )/ww
        DYYVs[k,l1,2] = sum(  (Y2s[I]                 - DYYs_bar[k,l1,2])^2 * taus[I,k] )/ww
        DYYVs[k,l1,3] = sum(  (Y3s[I] - B32s * Y2s[I] - DYYs_bar[k,l1,3])^2 * taus[I,k] )/ww
        DYYVs[k,l1,4] = sum(  (Y4s[I] - B43  * Y3s[I] - DYYs_bar[k,l1,4])^2 * taus[I,k] )/ww
      }
    }

     fitv  = slm.wfitc(XXV,c(as.numeric(DYYVm),as.numeric(DYYVs)),WW,CSw)$solution
     is=1
     S12   = sqrt(t(rdim(fitv[is:(is + nk*nf-1)],nk,nf))); is = is+nk*nf
     S2m   = sqrt(t(rdim(fitv[is:(is + nk*nf-1)],nk,nf))); is = is+nk*nf
     S3m   = sqrt(t(rdim(fitv[is:(is + nk*nf-1)],nk,nf))); is = is+nk*nf
     S43   = sqrt(t(rdim(fitv[is:(is + nk*nf-1)],nk,nf))); is = is+nk*nf
     S2s   = sqrt(t(rdim(fitv[is:(is + nk*nf-1)],nk,nf))); is = is+nk*nf
     S3s   = sqrt(t(rdim(fitv[is:(is + nk*nf-1)],nk,nf))); is = is+nk*nf
     if (ctrl$check_lik) {
        m4.mixt.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,A12=A12,A2ma=A2ma,A3ma=A3ma,A43=A43,A2mb=A2mb,A3mb=A3mb,S12=S12,S2m=S2m,S3m=S3m,S43=S43)
     }

    } else {
      # estimate the rhos -- use the other parameters
      # we start by recovering the posterior weight, and the variances for each term
      rwm  = c(t(taum + ctrl$posterior_reg))
      wv12 = c(t(S12[J1m,]))^2
      wv2  = c(t(S2m[J1m,]))^2
      wv32 = c(t(S3m[J2m,]))^2
      wv43 = c(t(S43[J2m,]))^2

      if (ctrl$est_rho[1]==TRUE) {
        fit = lm.wfit( as.matrix(c(t(Y2m - A2ma[J1m,]))), c(t(Y1m - A12[J1m,])),rwm/wv12)
        B12 = coefficients(fit)
        if (ctrl$check_lik) { m4.mixt.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,B12=B12)}
      }
      if (ctrl$est_rho[3]==TRUE) {
        fit = lm.wfit( as.matrix(c(t(Y3m - A3ma[J2m,]))), c(t(Y4m - A43[J2m,])),rwm/wv43)
        B43 = coefficients(fit)
        if (ctrl$check_lik) { m4.mixt.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,B43=B43)}
      }
      if (ctrl$est_rho[2]==TRUE) {
        fit = lm.wfit( as.matrix(c(t(Y2m - A2ma[J1m,] - A2mb[J2m]))), c(t(Y3m - A3ma[J2m,] - A3mb[J1m])),rwm/wv32)
        B32m = coefficients(fit)
        if (ctrl$check_lik) { m4.mixt.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,B32m=B32m)}
      }

      if (ctrl$est_rho[5]==TRUE) {
        rws  = c(t(taus + ctrl$posterior_reg))
        wv32 = c(t(S3s[J1s,]))^2
        fit     = lm.wfit( as.matrix(c(t(Y2s - A2s[J1s,]))), c(t(Y3s - A3s[J1s,])),rws/wv32)
        B32s    = coefficients(fit)
        if (ctrl$check_lik) m4.mixt.check.stayers(Y1s,Y2s,Y3s,Y4s,J1s,nk,Ns,model1,B32s=B32s);
      }
    }

    tic("mstep-ols")

    ## -------- PK probabilities ------------ #
    ## --- movers --- #
    for (l1 in 1:nf) for (l2 in 1:nf) {
      jj = l1 + nf*(l2-1)
      I = which(JJm == jj)
      if (length(I)>1) {
        pk1[jj,] = colSums(taum[I,])
      } else if (length(I)==0) { # this deals with the case where the cell is empty
        pk1[jj,] = 1/nk
      } else {
        pk1[jj,] = taum[I,]
      }
      pk1[jj,] = (pk1[jj,] + dprior-1 )/(sum(pk1[jj,] + dprior -1 ))
    }
    if (ctrl$check_lik) { m4.mixt.check(Y1m,Y2m,Y3m,Y4m,J1m,J2m,JJm,nk,Nm,model1,pk1=pk1)}

    ## --- stayers ----- #
    for (l1 in 1:nf) {
      I = which(J1s == l1)
      if (length(I)>1) {
        pk0[l1,] = colSums(taus[I,])
      } else if (length(I)==0) { # this deals with the case where the cell is empty
        pk0[l1,] = 1/nk
      } else {
        pk0[l1,] = taus[I,]
      }
      pk0[l1,] = (pk0[l1,] + dprior-1 )/(sum(pk0[l1,] + dprior -1 ))
    }

    #check_lik = computeLik(Y1m,Y2m,Y3m,Y4m,A12,B12,S12,A43,B43,S43,A2ma,A2mb,S2m,A3ma,A3mb,B32m,S3m)
    #if (check_lik<lik) cat("lik did not go down on pk1 update\n")

    tic("mstep-pks")

    # checking model fit
    if ((!any(is.na(model0))) & ((step %% ctrl$nplot) == (ctrl$nplot-1))) {
      I1 = order(colSums(A2ma))
      I2 = order(colSums(model0$A2ma))
      rr = addmom(A2ma[,I1],model0$A2ma[,I2],"A2ma")
      rr = addmom(A2mb,model0$A2mb,"A2mb",rr)
      rr = addmom(A3ma[,I1],model0$A3ma[,I2],"A3ma",rr)
      rr = addmom(A3mb,model0$A3mb,"A3mb",rr)
      rr = addmom(A2s[,I1],model0$A2s[,I2],"A2s",rr)
      rr = addmom(A12[,I1],model0$A12[,I2],"A12",rr)
      rr = addmom(A43[,I1],model0$A43[,I2],"A43",rr)
      rr = addmom(A3s[,I1],model0$A3s[,I2],"A3s",rr)
      rr = addmom(c(B12,B32m,B43), c(model0$B12,model0$B32m,model0$B43), "rhos", rr,type="rho")
      rr = addmom(S2m[,I1], model0$S2m[,I2], "S2m", rr,type="var")
      rr = addmom(S3m[,I1], model0$S3m[,I2], "S3m", rr,type="var")
      rr = addmom(S12[,I1],model0$S12[,I2],"S12",rr,type="var")
      rr = addmom(S43[,I1],model0$S43[,I2],"S43",rr,type="var")

      rr = addmom(pk1,model0$pk1,"pk1",rr,type="pr")
      rr = addmom(pk0,model0$pk0,"pk0",rr,type="pr")

      print(ggplot(rr,aes(x=val2,y=val1,color=type)) + geom_point() + facet_wrap(~name,scale="free") + theme_bw() + geom_abline(linetype=2))
    } else {
      if ((step %% ctrl$nplot) == (ctrl$nplot-1)) {
        mm = list(A12=A12,A2ma=A2ma,A3ma=A3ma,A43=A43)
        plot.endo.movers(mm)
      }
    }

    # -------- check convergence ------- #
    dlik = (lik_tot - lik_old)/abs(lik_old)
    lik_old = lik_tot
    lik_best = pmax(lik_best,lik_tot)
    if ( (step %% ctrl$ncat) == 0) flog.info("[%3i] levels=%i lik=%2.2f dlik=%2.2e dlik0=%2.2e liks=%2.2e likm=%2.2e r12=%.3f r23m=%.3f r23s=%.3f r43=%.3f",step,update_levels+0,lik_tot,dlik,lik_best-lik_tot,liks,likm,B12,B32m,B32s,B43);
    if (step>10) if (abs(dlik)<ctrl$tol) break;

    tic("loop-wrap")
  }

  # Y1 | Y2
  model$A12  = A12
  model$B12  = B12
  model$S12  = S12
  model$A43  = A43
  model$B43  = B43
  model$S43  = S43
  ## --movers --
  model$pk1  = pk1
  model$A2ma = A2ma
  model$A2mb = A2mb
  model$S2m  = S2m
  model$A3ma = A3ma
  model$A3mb = A3mb
  model$S3m  = S3m
  model$B32m = B32m
  # stayers
  model$pk0=pk0
  model$A2s=A2s
  model$A3s=A3s
  model$S2s=S2s
  model$S3s=S3s
  model$B32s=B32s

  model$NNm = acast(jdatae[,.N,list(j1,j2)],j1~j2,fill=0,value.var="N")
  model$lik=lik_tot
  model$lik_movers=likm
  model$lik_stayers=liks

  end.time <- Sys.time()
  time.taken <- end.time - start.time

  return(list(tic = tic(), model=model,lik=lik,step=step,dlik=dlik,time.taken=time.taken,ctrl=ctrl,liks=liks,likm=likm))
}

#' check the fit in the movers/stayers using imputed data
#' @export
m4.movers.checkfit <- function(jdata,r1,r4) {
  dd = jdata[, {
    d1=data.frame(src="data",
                  m1=mean(y1),m2=mean(y2),m3=mean(y3),m4=mean(y4),
                  d14=mean(y1-y4),m12=mean(y1-r1*y2_imp),m43=mean(y4-r4*y3),
                  cov12=cov(y1,y2),v1=var(y1),v2=var(y2))
    d2=data.frame(src="imp",
                  m1=mean(y1_imp),m2=mean(y2_imp),m3=mean(y3_imp),m4=mean(y4_imp),
                  d14=mean(y1_imp-y4_imp),m12=mean(y1_imp-r1*y2_imp),m43=mean(y4_imp-r4*y3_imp),
                  cov12=cov(y1_imp,y2_imp),v1=var(y1_imp),v2=var(y2_imp))
    rbind(d1,d2)
  },list(j1,j2)]

  ddm = melt(dd,id.vars = c("j1","j2","src"))
  ddm = cast(ddm,j1+j2+variable~src,value = "value")
  ggplot(ddm,aes(x=data,y=imp)) + geom_point() + facet_wrap(~variable,scales = "free") + theme_bw() +
     geom_abline(linetype=2)
  ddm
}

#' check the fit in the movers/stayers using imputed data
#' @export
m4.stayers.checkfit <- function(sdata,r1,r4) {
  dd = sdata[, {
    d1=data.frame(src="data",
                  m1=mean(y1),m2=mean(y2),m3=mean(y3),m4=mean(y4),
                  d14=mean(y1-y4),m12=mean(y1-r1*y2_imp),m43=mean(y4-r4*y3),
                  cov12=cov(y1,y2),v1=var(y1),v2=var(y2))
    d2=data.frame(src="imp",
                  m1=mean(y1_imp),m2=mean(y2_imp),m3=mean(y3_imp),m4=mean(y4_imp),
                  d14=mean(y1_imp-y4_imp),m12=mean(y1_imp-r1*y2_imp),m43=mean(y4_imp-r4*y3_imp),
                  cov12=cov(y1_imp,y2_imp),v1=var(y1_imp),v2=var(y2_imp))
    rbind(d1,d2)
  },list(j1)]

  ddm = melt(dd,id.vars = c("j1","src"))
  ddm = cast(ddm,j1+variable~src,value = "value")
  ggplot(ddm,aes(x=data,y=imp)) + geom_point() + facet_wrap(~variable,scales = "free") + theme_bw() +
    geom_abline(linetype=2)

  ddm
}

#' Approximate reclassification
#' @export
m4.mixt.estimate.reclassify <- function(sim,maxiter=20,split_movers=FALSE,ctrl,cl) {

  # store the original cluster values for comparaison
  sim$sdata[,jm1:=j1]
  sim$jdata[,jm1:=j1]
  sim$jdata[,jm2:=j2]

  # define classification and estimation samples
  sim$jdata[,splitdata := rank(runif(.N)) - 0.5 +0.1*(runif(1)-0.5) <= .N/2,f1]
  if (split_movers==TRUE) {
    sim$jdata[,split_classification := splitdata]
    sim$jdata[,split_estimation     := !splitdata]
  } else {
    sim$jdata[,split_classification := TRUE]
    sim$jdata[,split_estimation := TRUE]
  }

  # run the first estimation
  res_mixt_start = m4.mixt.estimate.all(list(sdata=sim$sdata,jdata=sim$jdata[split_estimation==TRUE]),nk=6,ctrl,cl)
  model = res_mixt_start$model

  nk = model$nk
  nf = model$nf

  rr_mixt = list()
  rr_mixt[[paste(0)]] = res_mixt_start
  for (iter in 1:maxiter) {

    # compute the posterior probability on stayers
    sdata.pos = sim$sdata[,{
      likm = rep(0,model$nf)
      # iterate on firm types
      for (ii in 1:.N) {
        ltau     = log(model$pk0)
        lnorm2 = lognormpdf(y2[ii]                                    , model$A2s, model$S2s)
        lnorm3 = lognormpdf(y3[ii] - model$B32s*(y2[ii] - model$A2s)  , model$A3s, model$S3s)
        lnorm1 = lognormpdf(y1[ii] - model$B12 *(y2[ii] - model$A2ma) , model$A12, model$S12)
        lnorm4 = lognormpdf(y4[ii] - model$B43 *(y3[ii] - model$A3ma) , model$A43, model$S43)
        lall = ltau + lnorm2 + lnorm4 + lnorm1  + lnorm3
        likm     = likm + blmrep:::logRowSumExp(lall)
      }
      # draw the type of the firm
      list(jp=1:model$nf,lik=likm,N=rep(.N,model$nf),jo=0)
    },list(fid=f1,jt=j1,jm=jm1)]

    # compute the posterior probability on movers
    jdata.pos.p1 = sim$jdata[split_classification==TRUE,{
      likm = rep(0,model$nf)
      # iterate on firm types
      ltau     = log(rdim(model$pk1,model$nf,model$nf,model$nk)[,jo,])
      A2mb     = array(model$A2mb[jo],c(nf,nk))
      A2ma     = model$A2ma
      A3ma     = spread(model$A3ma[jo,],1,nf)
      A3mb     = spread(model$A3mb,2,nk)
      A43      = spread(model$A43[jo,],1,nf)
      A12      = model$A12
      S2m      = model$S2m
      S12      = model$S12
      S3m      = spread(model$S3m[jo,],1,nf)
      S43      = spread(model$S43[jo,],1,nf)

      for (ii in 1:.N) {
        lnorm2   = lognormpdf(y2[ii]                                     , A2ma + A2mb     , S2m)
        lnorm1   = lognormpdf(y1[ii] - model$B12 *(y2[ii] - A2ma)        , A12             , S12)
        lnorm3   = lognormpdf(y3[ii] - model$B32m*(y2[ii] - A2ma - A2mb) , A3ma + A3mb     , S3m)
        lnorm4   = lognormpdf(y4[ii] - model$B43 *(y3[ii] - A3ma)        , A43             , S43)
        lall     = ltau + lnorm2 + lnorm4 + lnorm1  + lnorm3
        likm     = likm + blmrep:::logRowSumExp(lall)
      }
      # store everything
      list(jp=1:model$nf,lik=likm,N=rep(.N,model$nf))
    },list(fid=f1,jt=j1,jo=j2,jm=jm1)]

    jdata.pos.p2 = sim$jdata[split_classification==TRUE,{
      likm = rep(0,model$nf)
      # iterate on firm types
      ltau     = log(rdim(model$pk1,model$nf,model$nf,model$nk)[jo,,])
      A2mb     = spread(model$A2mb,2,nk)
      A2ma     = spread(model$A2ma[jo,],1,nf)
      A3ma     = model$A3ma
      A3mb     = array(model$A3mb[jo],c(nf,nk))
      A43      = model$A43
      A12      = spread(model$A12[jo,],1,nf)
      S2m      = spread(model$S2m[jo,],1,nf)
      S12      = spread(model$S12[jo,],1,nf)
      S3m      = model$S3m
      S43      = model$S43

      for (ii in 1:.N) {
        lnorm2   = lognormpdf(y2[ii]                                     , A2ma + A2mb     , S2m)
        lnorm1   = lognormpdf(y1[ii] - model$B12 *(y2[ii] - A2ma)        , A12             , S12)
        lnorm3   = lognormpdf(y3[ii] - model$B32m*(y2[ii] - A2ma - A2mb) , A3ma + A3mb     , S3m)
        lnorm4   = lognormpdf(y4[ii] - model$B43 *(y3[ii] - A3ma)        , A43             , S43)
        lall     = ltau + lnorm2 + lnorm4 + lnorm1  + lnorm3
        likm     = likm + blmrep:::logRowSumExp(lall)
      }
      # store everything
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
    res_mixt = m4.mixt.estimate.all(list(sdata=sim$sdata,jdata=sim$jdata[split_estimation==TRUE]),nk=6,ctrl,cl)

    model    = res_mixt$model
    rr_mixt[[paste(iter)]] = res_mixt
  }

  return(rr_mixt)
}




