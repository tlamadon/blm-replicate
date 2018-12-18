colSums2 <- function(M) {
  if (is.null(dim(M))) return(sum(M));
  if (length(dim(M))==1) return(sum(M));
  if (dim(M)[1]==1) return(MM[1,]);
  colSums(M)
}

m2.mixt.np.movers.lik <- function(Y1,Y2,J1,J2,model) {

  nk = model$nk
  Nm = length(Y1)
  nm = model$nm
  lpt      = array(0,c(Nm,nk,2,nm))  # log Pr[Y_it,m_it | k,t]
  post_m   = array(0,c(Nm,nk,2,nm))  # log Pr[ m   | k,t, data_i ])
  post_km  = array(0,c(Nm,nk,2,nm))  # log Pr[ m,k |   t, data_i ])
  post_k   = array(0,c(Nm,nk))       # log Pr[ k   |      data_i ])
  prior_k  = array(0,c(Nm,nk))       # prior probability
  lp       = array(0,c(Nm,nk))
  lk       = array(0,c(Nm,nk))
  lkt      = array(0,c(Nm,nk,2))
  liks = 0
  
  res = with(model,{
    liks = 0
    for (l1 in 1:nf) for (l2 in 1:nf) {
      I = which( (J1==l1) & (J2==l2))
      if (length(I)==0) next;
      
      for (k in 1:nk) {
        
        # For each (i,l,t) compute likelihood for each m|k
        for (m in 1:nm) {
          lpt[I,k,1,m] = lognormpdf(Y1[I], A1[l1,k] + C1[l1,k,m], S1[l1,k,m]) + log(W1[l1,k,m])  
          lpt[I,k,2,m] = lognormpdf(Y2[I], A2[l2,k] + C2[l2,k,m], S2[l2,k,m]) + log(W2[l2,k,m])  
        }
        
        # For each (i,l,k,t) compute likelihood, integrating m out
        # For each (i,l,k,t) compute likelihood, integrating m out
        for (t in 1:2) {
          lkt[I,k,t]     = logRowSumExp(lpt[I,k,t,])                                 
          post_m[I,k,t,] = exp(lpt[I,k,t,] - spread(logRowSumExp(lpt[I,k,t,]),2,nm)) # log(Pr[ m | k,t,Y_it ])
        }
        
        # For each (i) compute likelihood of k, summing periods
        lp[I,k]      = log(pk1[l1,l2,k]) + rowSums(lkt[I,k,])
        prior_k[I,k] = pk1[l1,l2,k]
      }
    }
    liks     = sum(logRowSumExp(lp))
    post_k   = exp(lp - spread(logRowSumExp(lp),2,nk)) # normalize the k probabilities Pr(k|Y1,Y2,Y3,Y4,l)
    list(post_k = post_k, post_m = post_m, lpt = lpt, lik= liks, prior_k=prior_k )
  }) 
  
  return(res)   
}

# computes Q and H for testing the EM algorithm
m2.mixt.np.movers.likhq <-  function(prior_k,lpt,post_k,post_m,post_k_bis,post_m_bis) {
  
  N  = dim(lpt)[1]
  nk = dim(post_k)[2]
  nm = dim(post_m)[4]
  
  # compute the likelihood
  llik = array(0,c(N,nk))
  for (k in 1:nk) {
    llik[,k] = log(prior_k[,k]) + logRowSumExp(lpt[,k,1,]) + logRowSumExp(lpt[,k,2,]) # logSumExp the mixtures, then sum the log-periods
  }  
  L = sum(logRowSumExp(llik))
  
  # compute Q
  Qtmp = array(0,c(N,nk,nm))
  for (k in 1:nk) {
    for (m in 1:nm) {
      Qtmp[,k,m] =  1/nm*post_k_bis[,k] * log(prior_k[,k]) +
                    post_k_bis[,k] * post_m_bis[,k,1,m] * lpt[,k,1,m] +
                    post_k_bis[,k] * post_m_bis[,k,2,m] * lpt[,k,2,m] 
    }
  }  
  Qv = sum(Qtmp)
  
  # compute H
  Htmp = array(0,c(N,nk,nm))
  for (k in 1:nk) {
    for (m in 1:nm) {
      Htmp[,k,m] =  1/nm*post_k_bis[,k] * log(post_k[,k]) +
                    post_k_bis[,k] * post_m_bis[,k,1,m] * log(post_m[,k,1,m]) +
                    post_k_bis[,k] * post_m_bis[,k,2,m] * log(post_m[,k,2,m])
    }
  }  
  Hv = -sum(Htmp)
  
  return(list(L=L,Q=Qv,H=Hv))
}

# here we want to check a bunch of properties for the EM steps
# model1 and model2 should be 2 consecutive steps
m2.mixt.np.movers.check <- function(Y1,Y2,J1,J2,model1,...) {
  
  change = list(...)

  model2 = copy(model1)
  model2[names(change)] = change[names(change)]

  res1 = m2.mixt.np.movers.lik(Y1,Y2,J1,J2,model1)
  res2 = m2.mixt.np.movers.lik(Y1,Y2,J1,J2,model2)

  res1$post_m[res1$post_m==0] = min(res1$post_m[res1$post_m>0])
  res2$post_m[res2$post_m==0] = min(res2$post_m[res2$post_m>0])
  
  LQH1 = m2.mixt.np.movers.likhq(res1$prior_k,res1$lpt,res1$post_k,res1$post_m,res1$post_k,res1$post_m)
  LQH2 = m2.mixt.np.movers.likhq(res2$prior_k,res2$lpt,res2$post_k,res2$post_m,res1$post_k,res1$post_m)
  
  # do the analysis, Evaluate Q(theta | theta^t) , Q(theta^t | theta^t), H(theta | theta^t) and H(theta^t | theta^t)
  # in both cases we integrate with respect to the posterior
  Q1 = LQH1$Q #sum( ( (res1$taum) * res1$lpm ))
  Q2 = LQH2$Q #sum( ( (res1$taum) * res2$lpm ))
  H1 = LQH1$H # - sum( (res1$taum) * log(res1$taum))
  H2 = LQH2$H #- sum( (res1$taum) * log(res2$taum))
  
  warn_str=""  
  test = TRUE
  if (is.finite(Q1*Q2*H1*H2)) {
    if (( Q2<Q1) | (H2<H1)) {
      warn_str = "!!!!!!!!!";
      test=FALSE
    }
  } else {
    flog.error("Na in Q or H")
  }
  flog.info("[emcheck] %s \tQd=%4.4e \tHd=%4.4e \tl0=%f \tl1=%f %s",paste(names(change),collapse = ","),  Q2-Q1,H2-H1,LQH1$L,LQH2$L,warn_str)
  return(test)
}



#' create a random model for EM with 
#' endogenous mobility with multinomial pr
#' @export
m2.mixt.np.new <-function(nk,nf,nm,inc=F,from.mini=NA,noise_scale=1) {
  
  model = list()

  # store the size
  model$nk    = nk
  model$nf    = nf
  model$nm    = nm

  # model for the mean of Y1,Y2
  model$A1    = array(0.9*(1 + 0.5*rnorm(nf*nk)),c(nf,nk))
  # model for Y4|Y3,l,k for movers and stayers
  model$A2    = array(0.9*(1 + 0.5*rnorm(nf*nk)),c(nf,nk))

  # model for p(K | l ,l') for movers
  model$pk1    = array(rdirichlet(nf*nf,rep(1,nk)),c(nf,nf,nk))
  # model for p(K | l ) for stayers
  model$pk0    = rdirichlet(nf,rep(1,nk))

  # model for the distribution  
  model$S1    = array(1+0.5*runif(nf*nk*nm),c(nf,nk,nm))*noise_scale   # variance of each component
  model$S2    = array(1+0.5*runif(nf*nk*nm),c(nf,nk,nm))*noise_scale
  model$C1    = array(rnorm(nm*nk*nf),c(nf,nk,nm))
  model$C2    = array(rnorm(nm*nk*nf),c(nf,nk,nm))
  model$W1    = array(rdirichlet(nf*nk,rep(1,nm)),c(nf,nk,nm))         # weight for each component
  model$W2    = array(rdirichlet(nf*nk,rep(1,nm)),c(nf,nk,nm))         
        
  # make sure the mixture is centered on 0
  for (l in 1:nf) for (k in 1:nk) {
    model$C1[l,k,] = model$C1[l,k,] - wtd.mean(model$C1[l,k,],model$W1[l,k,])
    model$C2[l,k,] = model$C2[l,k,] - wtd.mean(model$C2[l,k,],model$W2[l,k,])
  }

  # make the model increasing  
  if (inc) {
    for (l in 1:nf) {
      model$A1[l,]  = sort(model$A1[l,])
      model$A2[l,]  = sort(model$A2[l,])
    }
  }
  
  model$likm  =-Inf
  return(model)
}

#' create a random model for EM with 
#' endogenous mobility with multinomial pr
#' @export
m2.mixt.np.new.from.ns <-function(model,nm) {
  
  nf = model$nf
  nk = model$nk
  
  # add the number of mixtures
  model$nm    = nm
  model$pk1   = rdim(model$pk1,nf,nf,nk)
  
  # initiate the mixture components close to the normal model  
  model$S1    = spread(model$S1,3,nm)
  model$S2    = spread(model$S2,3,nm)
  model$C1    = spread(qnorm((1:nm)/(nm+1)),c(1,2),c(nf,nk))*mean(model$S1)
  model$C2    = spread(qnorm((1:nm)/(nm+1)),c(1,2),c(nf,nk))*mean(model$S2)
  model$W1    = array(1/nm,c(nf,nk,nm)) # set to equal and 1 half
  model$W2    = array(1/nm,c(nf,nk,nm))         
  
  return(model)
}



#' Simulate movers from the model
#' @export
m2.mixt.np.simulate.movers <- function(model,NNm) {
  
  J1 = array(0,sum(NNm)) 
  J2 = array(0,sum(NNm)) 
  Y1 = array(0,sum(NNm)) 
  Y2 = array(0,sum(NNm)) 

  M1 = array(0,sum(NNm)) 
  M2 = array(0,sum(NNm)) 

  K  = array(0,sum(NNm)) 
  
  A1  = model$A1
  A2  = model$A2

  S1 = model$S1
  S2 = model$S2
  C1 = model$C1
  C2 = model$C2
  W1 = model$W1
  W2 = model$W2

  pk1  = model$pk1
  nk   = model$nk
  nf   = model$nf
  nm   = model$nm
  
  i =1 
  for (l1 in 1:nf) for (l2 in 1:nf) {
    I = i:(i+NNm[l1,l2]-1)
    ni = length(I)
    jj = l1 + nf*(l2 -1)
    J1[I] = l1
    J2[I] = l2
    
    # draw k
    Ki = sample.int(nk,ni,T,pk1[l1,l2,])
    K[I] = Ki
    
    # draw the mixture components
    for (ik in 1:ni) {
        M1[I[ik]]     = sample.int(nm,1,prob=W1[l1,Ki[ik],],replace=TRUE)
        M2[I[ik]]     = sample.int(nm,1,prob=W2[l2,Ki[ik],],replace=TRUE)
        Y1[I[ik]]     = A1[l1,Ki[ik]] +   C1[l1,Ki[ik],M1[I[ik]]] +  S1[l1,Ki[ik],M1[I[ik]]] * rnorm(1)
        Y2[I[ik]]     = A2[l2,Ki[ik]] +   C2[l2,Ki[ik],M2[I[ik]]] +  S2[l2,Ki[ik],M2[I[ik]]] * rnorm(1)
    }

    i = i + NNm[l1,l2]
  }
  
  jdata = data.table(k=K,y1=Y1,y2=Y2,j1=J1,j2=J2,M1=M1,M2=M2)
  return(jdata)  
}

#' Simulate movers from the model
#' @export
m2.mixt.np.simulate.stayers <- function(model,NNs) {
  
  J1 = array(0,sum(NNs)) 
  J2 = array(0,sum(NNs)) 
  Y1 = array(0,sum(NNs)) 
  Y2 = array(0,sum(NNs)) 
  
  M1 = array(0,sum(NNs)) 
  M2 = array(0,sum(NNs)) 
  
  K  = array(0,sum(NNs)) 
  
  A1  = model$A1
  A2  = model$A2
  
  S1 = model$S1
  S2 = model$S2
  C1 = model$C1
  C2 = model$C2
  W1 = model$W1
  W2 = model$W2
  
  pk0  = model$pk0
  nk   = model$nk
  nf   = model$nf
  nm   = model$nm
  
  i =1 
  for (l1 in 1:nf) {
    I = i:(i+NNs[l1]-1)
    ni = length(I)
    J1[I] = l1
    J2[I] = l1
    
    # draw k
    Ki = sample.int(nk,ni,T,pk0[l1,])
    K[I] = Ki
    
    # draw the mixture components
    for (ik in 1:ni) {
      M1[I[ik]]     = sample.int(nm,1,prob=W1[l1,Ki[ik],],replace=TRUE)
      M2[I[ik]]     = sample.int(nm,1,prob=W2[l1,Ki[ik],],replace=TRUE)
      Y1[I[ik]]     = A1[l1,Ki[ik]] +   C1[l1,Ki[ik],M1[I[ik]]] +  S1[l1,Ki[ik],M1[I[ik]]] * rnorm(1)
      Y2[I[ik]]     = A2[l1,Ki[ik]] +   C2[l1,Ki[ik],M2[I[ik]]] +  S2[l1,Ki[ik],M2[I[ik]]] * rnorm(1)
    }
    
    i = i + NNs[l1]
  }
  
  sdata = data.table(k=K,y1=Y1,y2=Y2,j1=J1,j2=J2,M1=M1,M2=M2)
  return(sdata)  
}



# testing
m2.mixt.np.test <- function() {
  
  model = em.endo.level.np.new(3,4,3)
  Ns = rep(10000,4)
  sdata = em.endo.level.np.simulate.stayers.unc(model,Ns)
  
  # check that mean of the error is 0
  sdata[,list(mean(e1)/sd(e1),mean(e2)/sd(e2),mean(e3)/sd(e3),mean(e4)/sd(e4)),list(j1,k)]
  
  # check the mean of each component
  expect_true(sum( (acast(sdata[,mean(e2),list(j1,c2)],j1~c2,value.var = "V1")- model$C2s)^2)<0.01)
  expect_true(sum( (acast(sdata[,.N,list(j1,c2)][,p:=N/sum(N),j1],j1~c2,value.var = "p")- model$W2s)^2)<0.01)
  expect_true(sum( (acast(sdata[,mean(e1),list(j1,c1)],j1~c1,value.var = "V1")- model$C12)^2)<0.01)
  expect_true(sum( (acast(sdata[,mean(e3),list(j1,c3)],j1~c3,value.var = "V1")- model$C3s)^2)<0.01)
  expect_true(sum( (acast(sdata[,mean(e4),list(j1,c4)],j1~c4,value.var = "V1")- model$C43)^2)<0.01)
  expect_true(sum( (acast(sdata[,.N,list(j1,c1)][,p:=N/sum(N),j1],j1~c1,value.var = "p")- model$W12)^2)<0.01)
  expect_true(sum( (acast(sdata[,.N,list(j1,c3)][,p:=N/sum(N),j1],j1~c3,value.var = "p")- model$W3s)^2)<0.01)
  expect_true(sum( (acast(sdata[,.N,list(j1,c4)][,p:=N/sum(N),j1],j1~c4,value.var = "p")- model$W43)^2)<0.01)
  
  Ah = rdim(coef(sdata[,lm(y1 - model$B12*y2 ~ 0 + factor(k):factor(j1) )]),model$nf,model$nk)
  expect_true(  sum( ( model$A12 - model$B12 * model$A2s)^2)<0.1)
  
  expect_true(sum( (acast(sdata[,.N,list(j1,k)][,p:=N/sum(N),j1],j1~k,value.var = "p")- model$pk0)^2)<0.01)
  
}

#' Estimate a mixture of mixture model on stayers
#' @export
m2.mixt.np.stayers.estimate <- function(sdata,model,ctrl) {
  
  start.time <- Sys.time()
  tic <- tic.new()
  
  dprior = ctrl$dprior
  model0 = ctrl$model0
  tau    = ctrl$tau
  
  # ----- GET MODEL ---- #
  nk   = model$nk
  nf   = model$nf
  nm   = model$nm
  Nm   = sdata[,.N]

  A1  = model$A1
  A2  = model$A2
  
  S1 = model$S1
  S2 = model$S2
  C1 = model$C1
  C2 = model$C2
  W1 = model$W1
  W2 = model$W2

  pk0 = model$pk0
    
  # ----- GET DATA ---- #
  # movers
  Y1 = sdata$y1
  J1 = sdata$j1

  lik_old  = -Inf
  lik      = -Inf
  lik_best = -Inf
  likm=0
  dlik_ma = 1
  
  model$dprior = dprior
  model1 = copy(model)
  
  tic("prep")
  
  lpt     = array(0,c(Nm,nk,2,nm))  # log(Pr_i[k,Y1,t,m])
  post_m  = array(0,c(Nm,nk,2,nm))  # log(Pr_i[k,Y1,t,m])
  lp      = array(0,c(Nm,nk))
  lk      = array(0,c(Nm,nk))
  post_k  = array(0,c(Nm,nk))
  lkt     = array(0,c(Nm,nk,2))
  
  stop = F;
  for (step in 1:ctrl$maxiter) {
    
    model1[['pk0']] = pk0

    # ---------- E STEP ------------- #
    # compute the tau probabilities and the likelihood
    if (is.na(tau[1]) | (step>1)) {
      
      liks = 0
      for (l1 in 1:nf){
        I = which((J1==l1))
        if (length(I)==0) next;
        for (k in 1:nk) {
          for (m in 1:nm) {
            lpt[I,k,1,m] = lognormpdf(Y1[I], A1[l1,k] + C1[l1,k,m], S1[l1,k,m]) + log(W1[l1,k,m])  
          }
          # take the expectation on the mixture realization within each period
          for (t in 1) {
            lkt[I,k,t]     = logRowSumExp(lpt[I,k,t,])
            post_m[I,k,t,] = exp(lpt[I,k,t,] - spread(logRowSumExp(lpt[I,k,t,]),2,nm))
          }
          
          # sum the log of the periods 
          lp[I,k] = log(pk0[l1,k]) + lkt[I,k,1]
        }
      }
      liks     = sum(logRowSumExp(lp))
      post_k   = exp(lp - spread(logRowSumExp(lp),2,nk)) # normalize the k probabilities Pr(k|Y1,Y2,Y3,Y4,l)
      
      # compute prior
      lik_prior = (dprior-1) * sum(log(pk0))
      lik = liks + lik_prior
      
    } else {
      cat("skiping first max step, using supplied posterior probabilities\n")
    }
    
    tic("estep")
    
    if (stop) break;
    
    # ----- M-STEP ------- #
    
    # ====================== update the pk0 ===================== #
    for (l1 in 1:nf)  {
      I = which((J1==l1))
      if (length(I)>1) {
        pk0[l1,] = colSums(post_k[I,])
      } else {
        pk0[l1,] = post_k[I,]
      }
      pk0[l1,] = (pk0[l1,] + dprior-1)/sum(pk0[l1,]+dprior-1)
    }
    tic("mstep-pks")
    
    # checking model fit
    if ((!any(is.na(model0))) & ((step %% ctrl$nplot) == (ctrl$nplot-1))) {
      rr = addmom(A1,  model0$A1,  "A1")
      rr = addmom(A2,  model0$A2,  "A2", rr)
      rr = addmom(C1,  model0$C1,  "C1", rr)
      rr = addmom(C2,  model0$C2,  "C2", rr)
      rr = addmom(W1,  model0$W1,  "W1", rr)
      rr = addmom(W2,  model0$W2,  "W2", rr)
      rr = addmom(S1,  model0$S1,  "S1",  rr,type="var")
      rr = addmom(S2,  model0$S2,  "S2",  rr,type="var")
      rr = addmom(pk1, model0$pk1, "pk1", rr,type="pr")
      rr = addmom(pk0, model0$pk0, "pk0", rr,type="pr")
      print(ggplot(rr,aes(x=val1,y=val2,color=type)) + geom_point() + facet_wrap(~name,scale="free") + theme_bw() + geom_abline(linetype=2))
    } 
    
    # -------- check convergence ------- #
    dlik = (lik - lik_old)/abs(lik_old)
    lik_old = lik
    lik_best = pmax(lik_best,lik)
    dlik_smooth = 0
    flog.info("[%3i] lik=%4.4f dlik=%4.4e (%4.4e) dlik0=%4.4e liks=%4.4e likm=%4.4e",step,lik,dlik,dlik_smooth,lik_best-lik,liks,model$likm,name="m2.mixt.np")
    if (step>10) {
      dlik_smooth = 0.5*dlik_smooth + 0.5*dlik
      if (abs(dlik_smooth)<ctrl$tol) break;
    }
    
    tic("loop-wrap")
  }
  
  # -------- store results --------- #
  model$liks = liks
  model$pk0  = pk0
  model$NNs  = sdata[,.N,j1][order(j1)][,N]
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  return(list(tic = tic(), model=model,lik=lik,step=step,dlik=dlik,time.taken=time.taken,ctrl=ctrl,liks=liks,likm=likm))  
}

#' Estimate a mixture of mixture model 
#' @export
m2.mixt.np.movers.estimate <- function(jdata,model,ctrl) {
  
  start.time <- Sys.time()
  tic <- tic.new()
  
  dprior = ctrl$dprior
  model0 = ctrl$model0
  tau    = ctrl$tau
  
  # ----- GET MODEL ---- #
  nk   = model$nk
  nf   = model$nf
  nm   = model$nm
  nt   = 2
  Nm   = nrow(jdata)
  
  A1  = model$A1
  A2  = model$A2
  
  S1 = model$S1
  S2 = model$S2
  C1 = model$C1
  C2 = model$C2
  W1 = model$W1
  W2 = model$W2
  
  pk1 = model$pk1
  
  # ----- GET DATA ---- #
  # movers
  Y1 = jdata$y1
  Y2 = jdata$y2
  J1 = jdata$j1
  J2 = jdata$j2
  
  # constructing regression matrix
  XXA  = diag(2) %x% diag(nf) %x% diag(nk)  %x% rep(1,nm)  # this is for the A1,A2 coefficients
  XXC  = diag(2) %x% diag(nf) %x% diag(nk)  %x% diag(nm)   # this is for the C1,C2 coefficients
  XX   = cBind(XXA,XXC) # we combine them
  
  # prepare the constraints that force the mean of the within mixtures to be 0
  # for each l,k it has to sum to 0, the weights will be added in later
  CM  = cons.pad(cons.sum(nm,2*nf*nk),2*nk*nf,0) # we constraint the sum, we pad the begining because of the As
  
  # add the constraints from the user for the As
  CS1  = cons.pad(cons.get(ctrl$cstr_type[1],ctrl$cstr_val[1],nk,nf),nk*nf*0, nk*nf*(1+nt*nm))
  CS2  = cons.pad(cons.get(ctrl$cstr_type[2],ctrl$cstr_val[2],nk,nf),nk*nf*1, nk*nf*nt*nm)
  # combine them
  CS = cons.bind(CS1,CS2)
  
  # create the stationary contraints
  if (ctrl$fixb==T) {
    CS2 = cons.pad( cons.fixb(nk,nf,2),0, nk*nf*nt*nm)
    CS  = cons.bind(CS2,CS)
  }
  
  # create a constraint for the variances, here we have nm*nf variances
  if (ctrl$model_var==T) {
    CSw = cons.none(nm,nf*2)
  } else{
    CS1  = cons.pad(cons.mono_k(nm,nf),nm*nf*0, nm*nf*3)
    CS2  = cons.pad(cons.mono_k(nm,nf),nm*nf*1, nm*nf*2)
    CSw  = cons.bind(CS1,CS2)
    CSw$meq = length(CSw$H)
  }
  
  lik_old  = -Inf
  lik      = -Inf
  lik_best = -Inf
  likm=0
  dlik_ma = 1
  
  model$dprior = dprior
  model1 = copy(model)
  
  tic("prep")
  
  lpt     = array(0,c(Nm,nk,2,nm))  # log(Pr_i[k,Y1,t,m])
  post_m  = array(0,c(Nm,nk,2,nm))  # log(Pr_i[k,Y1,t,m])
  lp      = array(0,c(Nm,nk))
  lk      = array(0,c(Nm,nk))
  post_k  = array(0,c(Nm,nk))
  lkt     = array(0,c(Nm,nk,2))
  
  stop = F;
  for (step in 1:ctrl$maxiter) {
    
    model1[c('A1','A2','pk1','C1','C2','W1','W2','S1','S2')] = list(A1,A2,pk1,C1,C2,W1,W2,S1,S2)
    
    # ---------- E STEP ------------- #
    # compute the tau probabilities and the likelihood
    if (is.na(tau[1]) | (step>1)) {
      
      liks = 0
      for (l1 in 1:nf) for (l2 in 1:nf) {
        I = which( (J1==l1) & (J2==l2))
        if (length(I)==0) next;
        
        for (k in 1:nk) {
          for (m in 1:nm) {
            lpt[I,k,1,m] = lognormpdf(Y1[I], A1[l1,k] + C1[l1,k,m], S1[l1,k,m]) + log(W1[l1,k,m])  
            lpt[I,k,2,m] = lognormpdf(Y2[I], A2[l2,k] + C2[l2,k,m], S2[l2,k,m]) + log(W2[l2,k,m])  
          }
          # take the expectation on the mixture realization within each period
          for (t in 1:2) {
            lkt[I,k,t]     = logRowSumExp(lpt[I,k,t,])
            post_m[I,k,t,] = exp(lpt[I,k,t,] - spread(logRowSumExp(lpt[I,k,t,]),2,nm))
          }
          
          # sum the log of the periods 
          lp[I,k] = log(pk1[l1,l2,k]) + rowSums(lkt[I,k,])
        }
      }
      liks     = sum(logRowSumExp(lp))
      post_k   = exp(lp - spread(logRowSumExp(lp),2,nk)) # normalize the k probabilities Pr(k|Y1,Y2,Y3,Y4,l)
      
      # compute prior
      lik_prior = (dprior-1) * sum(log(pk1))
      lik = liks + lik_prior
      
    } else {
      cat("skiping first max step, using supplied posterior probabilities\n")
    }
    
    tic("estep")
    
    if (stop) break;
    
    # ----- M-STEP ------- #
    
    # ====================== update the pk1 ===================== #
    for (l1 in 1:nf) for (l2 in 1:nf) {
      I = which( (J1==l1) & (J2==l2))
      if (length(I)>1) {
        pk1[l1,l2,] = colSums(post_k[I,])
      } else {
        pk1[l1,l2,] = post_k[I,]
      }
      pk1[l1,l2,] = (pk1[l1,l2,] + dprior-1)/sum(pk1[l1,l2,]+dprior-1)
    }
    if (ctrl$check_lik) m2.mixt.np.movers.check(Y1,Y2,J1,J2,model1,pk1=pk1);
    
    update_levels = ((step%%2)==0) | (step<3) 
    
    if (FALSE) {
      
      
      
    } else {
      
      # ================== Mixture weights ==================== #
      # how to do this under constraints?
      for (l in 1:nf) {
        I1 = which(J1==l)
        I2 = which(J2==l)
        for (k in 1:nk) {
          W1[l,k,] = colSums2(post_m[I1,k,1,]*post_k[I1,k])
          W2[l,k,] = colSums2(post_m[I2,k,2,]*post_k[I2,k])
          W1[l,k,]  = (W1[l,k,] + dprior-1)/sum(W1[l,k,]+dprior-1)
          W2[l,k,]  = (W2[l,k,] + dprior-1)/sum(W2[l,k,]+dprior-1)
        }
      }
      if (ctrl$check_lik) m2.mixt.np.movers.check(Y1,Y2,J1,J2,model1,pk1=pk1,W1=W1,W2=W2);

      # ================== update means by constructing reduced matrix ==================== #
      for (l in 1:nf) {
        I1 = which(J1==l)
        I2 = which(J2==l)
        for (k in 1:nk) {
          # update the constraints with the new weights
          # this is used to estimate the means (the mixture are mean
          # 0 in the way we estimate them)
          # the rows run T x L X K
          # the cols run T x L x K x M ( and they start with 2*n*k A values)
          CM$C[ k + nk*( l-1 + nf*0) ,2*nf*nk + 1:nm + nm*(k-1 + nk*(l-1 + nf*0))] = W1[l,k,]
          CM$C[ k + nk*( l-1 + nf*1) ,2*nf*nk + 1:nm + nm*(k-1 + nk*(l-1 + nf*1))] = W2[l,k,]
        }
      }
      
      # this nt*nk*nl*nm equations and as many dummies, with the proper constraints
      # what we need is the weight on each of this equations
      ww = array(0,c(nm,nk,nf,nt))
      YY = array(0,c(nm,nk,nf,nt))
      for (l in 1:nf) {
        I1 = which( (J1==l))
        I2 = which( (J2==l))
        for (m in 1:nm) {
          for (k in 1:nk) {
            # T=1
            w = (post_k[I1,k] * post_m[I1,k,1,m]+ ctrl$posterior_reg)/S1[l,k,m]^2 
            ww[m,k,l,1] = sum(w)
            #YY[m,k,l,1] = sum((Y1[I1]-C1[l,k,m])*w)/ww[m,k,l,1]
            YY[m,k,l,1] = sum((Y1[I1])*w)/ww[m,k,l,1]
            
            # T=2
            w = (post_k[I2,k] * post_m[I2,k,2,m]+ ctrl$posterior_reg)/S2[l,k,m]^2 
            ww[m,k,l,2] = sum(w)
            #YY[m,k,l,2] = sum((Y2[I2] - C2[l,k,m]) * w)/ww[m,k,l,2]
            YY[m,k,l,2] = sum((Y2[I2]) * w)/ww[m,k,l,2]
          }
        }
      }
      
      # # checking YY
      # for (l in 1:nf) for (k in 1:nk) {
      #   cat(sprintf("[k=%i,l=%i] %f/%f\n",k,l,YY[1:nm + nm*(k-1 + nk*(l-1))], model$A1[l,k] + model$C1[l,k,]))
      # }
      
      CS2 = cons.bind(CS,CM)
      #CS2 = CS
      #CS2$C = CS2$C[,1:ncol(XXA)]
      # then we regress
      # flog.info("min S1=%e S2=%e s=%f",min(S1),min(S2),min(ww))
      fit  = slm.wfitc(XX,as.numeric(YY),as.numeric(ww),CS2)$solution
      A1   = t(rdim(fit[1:(nk*nf)],nk,nf))     
      A2   = t(rdim(fit[(nk*nf+1):(2*nk*nf)],nk,nf))
      for (l in 1:nf) {
        for (k in 1:nk) {
          C1[l,k,]   = fit[ 2*nk*nf + 1:nm + nm*(k-1 + nk*(l-1 + nf*0))]
          C2[l,k,]   = fit[ 2*nk*nf + 1:nm + nm*(k-1 + nk*(l-1 + nf*1))]
        }
      }
      pred = rdim(XX %*% fit,nm,nk,nf,2)
      if (ctrl$check_lik) m2.mixt.np.movers.check(Y1,Y2,J1,J2,model1,pk1=pk1,A1=A1,A2=A2,W1=W1,W2=W2,C1=C1,C2=C2);
      
      # update variances by contructing aggregated residuals    
      ww = array(0,c(nm,nk,nf,nt))
      YY = array(0,c(nm,nk,nf,nt))
      for (l in 1:nf) {
        I1 = which( (J1==l))
        I2 = which( (J2==l))
        for (m in 1:nm) {
          for (k in 1:nk) {
            # T=1
            w = (post_k[I1,k] * post_m[I1,k,1,m] + ctrl$posterior_reg)/S1[l,k,m]^2
            ww[m,k,l,1] = sum(w)
            YY[m,k,l,1] = sum((Y1[I1]-pred[m,k,l,1])^2*w)/ww[m,k,l,1]
            
            # T=2
            w = (post_k[I2,k] * post_m[I2,k,2,m] + ctrl$posterior_reg)/S2[l,k,m]^2
            ww[m,k,l,2] = sum(w)
            YY[m,k,l,2] = sum((Y2[I2]-pred[m,k,l,2])^2 * w)/ww[m,k,l,2]
          }
        }
      }
      
      
      fitv  = lm.wfitnn(XXC,as.numeric(YY),as.numeric(ww))$solution
      for (l in 1:nf) {
        for (k in 1:nk) {
          S1[l,k,]   = sqrt(fitv[ 1:nm + nm*(k-1 + nk*(l-1 + nf*0))])
          S2[l,k,]   = sqrt(fitv[ 1:nm + nm*(k-1 + nk*(l-1 + nf*1))])
        }
      }
      if (ctrl$check_lik) m2.mixt.np.movers.check(Y1,Y2,J1,J2,model1,pk1=pk1,W1=W1,W2=W2,S1=S1,S2=S2,A1=A1,A2=A2,C1=C1,C2=C2);
      S1[S1<ctrl$sd_floor]=ctrl$sd_floor # having a variance of exacvtly 0 creates problem in the likelihood
      S2[S2<ctrl$sd_floor]=ctrl$sd_floor
    }
    
    tic("mstep-pks")
    
    # checking model fit
    if ((!any(is.na(model0))) & ((step %% ctrl$nplot) == (ctrl$nplot-1))) {
      rr = addmom(A1,  model0$A1,  "A1")
      rr = addmom(A2,  model0$A2,  "A2", rr)
      rr = addmom(C1,  model0$C1,  "C1", rr)
      rr = addmom(C2,  model0$C2,  "C2", rr)
      rr = addmom(W1,  model0$W1,  "W1", rr)
      rr = addmom(W2,  model0$W2,  "W2", rr)
      rr = addmom(S1,  model0$S1,  "S1",  rr,type="var")
      rr = addmom(S2,  model0$S2,  "S2",  rr,type="var")
      rr = addmom(pk1, model0$pk1, "pk1", rr,type="pr")
      print(ggplot(rr,aes(x=val1,y=val2,color=type)) + geom_point() + facet_wrap(~name,scale="free") + theme_bw() + geom_abline(linetype=2))
    } 
    
    # -------- check convergence ------- #
    dlik = (lik - lik_old)/abs(lik_old)
    lik_old = lik
    lik_best = pmax(lik_best,lik)
    dlik_smooth = 0
    flog.info("[%3i] lik=%4.4f dlik=%4.4e (%4.4e) dlik0=%4.4e liks=%4.4e likm=%4.4e",step,lik,dlik,dlik_smooth,lik_best-lik,liks,model$likm,name="m2.mixt.np")
    if (step>10) {
      dlik_smooth = 0.5*dlik_smooth + 0.5*dlik
      if (abs(dlik_smooth)<ctrl$tol) break;
    }
    
    tic("loop-wrap")
  }
  
  # -------- store results --------- #
  model$liks = liks
  model[c('A1','A2','S1','S2','pk1','C1','C2','W1','W2')] = list(A1,A2,S1,S2,pk1,C1,C2,W1,W2)

  model$NNm  = acast(jdata[,.N,list(j1,j2)],j1~j2,fill=0,value.var="N")
  model$likm = liks  
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  return(list(tic = tic(), model=model,lik=lik,step=step,tau=tau,dlik=dlik,time.taken=time.taken,ctrl=ctrl,liks=liks,likm=likm))  
}


#' Compute some moments on the distribution of residuals
#'
#' export
m2.mixt.np.movers.residuals <- function(model,plot=T,n=100) {
  W1= model$W1
  S1= model$S1
  C1= model$C1
  
  nk= model$nk
  nf= model$nf
  nm= model$nm
  
  m1 = array(0,c(nk,nf))
  m2 = array(0,c(nk,nf))
  m3 = array(0,c(nk,nf))
  m4 = array(0,c(nk,nf))
  
  for (l in 1:nf) 
    for(k in 1:nk) {
      # compute the mena/variance/curtosis/skewness of the error distributions
      m1[k,l] = sum(W1[l,k,]*C1[l,k,])
      m2[k,l] = sqrt(sum( W1[l,k,]*(C1[l,k,]^2+S1[l,k,]^2)) - m1[k,l])
      m3[k,l] = 1/m2[k,l]^3 * sum( W1[l,k,]*(C1[l,k,] - m1[k,l])*(3*S1[l,k,]^2 +  ( C1[l,k,] - m1[k,l]   )^2))
      m4[k,l] = 1/m2[k,l]^4 * sum( W1[l,k,]*( 3*S1[l,k,]^4 +  6*( C1[l,k,] - m1[k,l]   )^2*S1[l,k,]^2 + ( C1[l,k,] - m1[k,l]   )^4))
    }
  
  rr = data.frame()
  for (l in 1:nf) 
    for(k in 1:nk) {
      
      # get the standard deviation
      msd = sqrt(sum( W1[l,k,]*(C1[l,k,]^2+S1[l,k,]^2)))
      # cover 99% 
      xsup = qnorm((1:n)/(n+1),sd=msd)
      # store values    
      ysup = 0
      for (m in 1:nm) {
        ysup = ysup +  W1[l,k,m] * dnorm( xsup, mean = C1[l,k,m], sd = S1[l,k,m]) 
      }
      ysup = ysup/max(ysup)
      rr = rbind(rr,data.frame(l=l,k=k,x= xsup,y=ysup))
      
      # compute the mena/variance/curtosis/skewness of the error distribution
      m1[k,l] = 0
      m2[k,l] = msd
      m3[k,l] = 1/m2[k,l]^3 * sum( W1[l,k,]*(C1[l,k,] - m1[k,l])*(3*S1[l,k,]^2 +  ( C1[l,k,] - m1[k,l]   )^2))
      m4[k,l] = 1/m2[k,l]^4 * sum( W1[l,k,]*( 3*S1[l,k,]^4 +  6*( C1[l,k,] - m1[k,l]   )^2*S1[l,k,]^2 + ( C1[l,k,] - m1[k,l]   )^4))
  }


  if (plot==TRUE) {
    print(ggplot(rr,aes(x=x,y=y)) + geom_line() + facet_grid(k~l,scales = "free") + theme_bw())
  }
  
  return(list(mean=m1,sd=m2,skew=m3,kurt=m4,all=rr))
}



test.em.level.np <- function() {
  rm(list = ls())

  # testing likelihood
  model = m2.mixt.np.new(nk = 2,nf = 4,nm = 3,noise_scale = 1,inc=T)
  Nm    = round(array(1000/(model$nf*model$nk),c(model$nf,model$nf)))
  jdata = m2.mixt.np.simulate.movers(model,Nm)
  res   = blm:::m2.mixt.np.movers.lik(jdata$y1,jdata$y2,jdata$j1,jdata$j2,model)
  LQH1  = blm:::m2.mixt.np.movers.likhq(res$prior_k,res$lpt,res$post_k,res$post_m,res$post_k,res$post_m)
  (LQH1$Q + LQH1$H)/LQH1$L
  
  
  # Estimation
  model = m2.mixt.np.new(nk = 2,nf = 4,nm = 3,noise_scale = 0.1,inc=T)
  model_rand = m2.mixt.np.new(nk = model$nk,nf = model$nf,nm = model$nm,noise_scale = 0.3,inc=T)
  Nm = round(array(20000/(model$nf*model$nk),c(model$nf,model$nf)))
  jdata = m2.mixt.np.simulate.movers(model,Nm)

  model_bis = copy(model)
  model_bis$A1=model_rand$A1
  model_bis$A2=model_rand$A2
  model_bis$pk1=model_rand$pk1
  model_bis$W1=model_rand$W1
  model_bis$W2=model_rand$W2
  
  ctrl   = em.control(nplot=10,check_lik=T,fixb=F,est_rho=F,model0=model,dprior=1.05)
  res    = m2.mixt.np.movers.estimate(jdata,model_bis,ctrl)
  
  rr = blm:::m2.mixt.np.movers.residuals(model)
  rr = blm:::m2.mixt.np.movers.residuals(model,TRUE)
   
}


#' Computes the variance decomposition by simulation
#' @export
m2.mixt.np.vdec <- function(model,nsim,stayer_share=1,ydep="y1") {
  
  if (ydep!="y1") flog.warn("ydep other than y1 is not implemented, using y1")
  # simulate movers/stayers, and combine
  NNm = model$NNm
  NNs = model$NNs
  NNm[!is.finite(NNm)]=0
  NNs[!is.finite(NNs)]=0
  
  NNs = round(NNs*nsim*stayer_share/sum(NNs))
  NNm = round(NNm*nsim*(1-stayer_share)/sum(NNm))
  flog.info("computing var decomposition with ns=%i nm=%i",sum(NNs),sum(NNm))

  # we simulate from the model both movers and stayers
  sdata.sim = m2.mixt.np.simulate.stayers(model,NNs)
  jdata.sim = m2.mixt.np.simulate.movers(model,NNm)
  sdata.sim = rbind(sdata.sim[,list(j1,k,y1)],jdata.sim[,list(j1,k,y1)])
  proj_unc  = lin.proj(sdata.sim,"y1","k","j1");  
  
  return(proj_unc)
}
 