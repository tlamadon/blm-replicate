# This implements the probabilistic model (in stationary)

# ------------- Initiliazing functions ---------------------

# create a random model for EM with
# endogenous mobility with multinomial pr
#' @export
m2.proba.new <-function(nk,nf,fixb=F) {

  model = list()
  # model for Y1,Y2
  model$A1    = array(0.9*(1 + 0.5*rnorm(nf*nk)),c(nf,nk))
  model$S1    = 0.3*array(1+0.5*runif(nf*nk),c(nf,nk))
  model$A2    = array(0.9*(1 + 0.5*rnorm(nf*nk)),c(nf,nk))
  model$S2    = 0.3*array(1+0.5*runif(nf*nk),c(nf,nk))

  # model for p(K | l ,l') for movers
  model$g1    = rdim(rdirichlet(nf*nf,rep(1,nk)),nf,nf,nk)
  # model for p(K | l ,l') for stayers
  model$g0    = rdirichlet(nf,rep(1,nk))
  # model for mobility
  model$pm    = array(0.1,c(nf,nk))
  # distribution for k and l
  model$pk    = rdirichlet(1,rep(1,nk))
  model$pl    = rdirichlet(1,rep(1,nf))

  model$nk    = nk
  model$nf    = nf

  for (l in 1:nf) {
   model$A1[l,] = sort(model$A1[l,])
  }

  return(model)
}

#' create proba model from mixt2 model
#'
#' @export
m2.proba.new.from.mixt <-function(model0,sdata_subsample=0.1,mobility_smooth=0.001) {

  model = list()
  # wages are identical
  model[c('A1','A2','S1','S2','nk','nf')]    = model0[c('A1','A2','S1','S2','nk','nf')]

  nf = model0$nf
  nk = model0$nk

  pk0 = model0$pk0[1,,]
  pk1 = rdim(model0$pk1,model0$nf,model0$nf,model0$nk)
  Ns  = model0$NNs/sdata_subsample + as.numeric(rowSums(model0$NNm))

  pr_joint = (pk0+mobility_smooth) * spread(Ns/sum(Ns),2,model0$nk)

  # distribution for k
  model$pk    = colSums(pr_joint)

  # model for choice of class in period 1
  model$g0    = log(pr_joint / spread(model$pk,1,model0$nf))

  # model for mobility
  model$pm = array(0,c(nf,nk))
  for (l in 1:nf) for (k in 1:nk) {
    model$pm[l,k] = sum(model0$NNm[l,] * pk1[l,,k])/(model0$NNs[l]/sdata_subsample * pk0[l,k] + sum(model0$NNm[l,] * pk1[l,,k]))
  }

  # model for choice of class in period 2 conditional on moving
  model$g1 = array(0,c(nf,nf,nk))
  pr_joint1 = (pk1) * spread((model0$NNm+mobility_smooth)/sum(model0$NNm+mobility_smooth),3,model0$nk)
  for (l1 in 1:nf) for (k in 1:k) {
    model$g1[l1,,k] = log(pr_joint1[l1,,k]/sum(pr_joint1[l1,,k]))
  }

  return(model)
}





# ------------- Simulating functions ---------------------

#' impute wages on existing data
#' @export
m2.proba.impute.stayers <- function(sdatae,model) {

  A1  = model$A1
  S1  = model$S1
  g0  = model$g0
  A2  = model$A2
  S2  = model$S2
  nk  = model$nk
  nf  = model$nf

  # we compute the probability of k|j1
  pk0 = exp(model$g0) * spread(as.numeric(model$pk),1,model$nf)
  pk0 = pk0 / spread(rowSums(pk0),2,model$nk)

  # ------  impute K, Y1, Y4 on jdata ------- #
  sdatae.sim = copy(sdatae)
  sdatae.sim[, c('k_imp','y1_imp') := {
    ni = .N
    Ki  = sample.int(nk,.N,prob = pk0[j1,],replace=T)
    # draw Y2, Y3
    Y1  = A1[j1,Ki] + S1[j1,Ki] * rnorm(ni)
    list(Ki,Y1)
  },list(j1)]

  return(sdatae.sim)
}

# -------------------- Estimating functions -----------------------------

#' computing the likelihood of stayers for a given group
#' of workers in a cluster
#'
#' @param Y1 vector of log-wages in period 1
#' @param J1s vector of firm classes in period 1
#' @param sj1 a vector a log-size of the firm when the worker is employed
#' @param model a probabilistic model
#' @param lik.den pre-computed denominators
#' @export
m2.proba.lik.stayers <- function(Y1,J1s,sj1,model,lik.den,include_j=T) {
  lps = array(0,c(length(Y1),model$nk))

  for (l1 in 1:model$nf)  {
    I = which(J1s==l1)
    for (k in 1:model$nk) {
      lpr_y1  = lognormpdf(Y1[I] , model$A1[l1,k], model$S1[l1,k])  # the wage
      lpr_l1  = model$g0[l1,k] - lik.den$R0[k]                      # choosing the current cluster given the type k
      lpr_j1  = sj1[I] - lik.den$Rj[l1]                                 # choosing the current firm within the cluster
      lps[I,k] = log(model$pk[k]) + lpr_l1+ lpr_j1 + lpr_y1 + log(1-model$pm[l1,k])
    }
  }

  sum(logRowSumExp(lps))
}

#' computing the likelihood of stayers for a given group
#' of workers in a cluster
#'
#' @param Y1 vector of log-wages in period 1
#' @param Y2 vector of log-wages in period 2
#' @param J1m vector of firm types in period 1
#' @param J2m vector of firm types in period 2
#' @param sj1 a vector a log-size of the firm when the worker is employed in period 1
#' @param sj2 a vector a log-size of the firm when the worker is employed in period 2
#' @param R0 a vector of constants in the probabilistic model for choice of cluster in period 1
#' @param R1 a vector of constants in the probabilistic model for choice of cluster in period 2
#' @param Rj1 a vector of constants in the probabilistic model for choice of a firm within a cluster in period 1
#' @param Rj2 a vector of constants in the probabilistic model for choice of a firm within a cluster in period 2
#' @param model a probabilistic model
#' @export
m2.proba.lik.movers <- function(Y1m,Y2m,J1m,J2m,sj1,sj2,model,lik.den,include_j=T) {
  lpm = array(0,c(length(Y1m),model$nk))

  for (l1 in 1:model$nf) for (l2 in 1:model$nf) {
    I = which( (J1m==l1) & (J2m==l2))
    if (length(I)==0) next;

    for (k in 1:model$nk) {
      # worker mobility
      lpr_k1    = model$g0[l1,k] - lik.den$R0[k]
      lpr_move  = log(model$pm[l1,k])
      lpr_k2    = model$g1[l1,l2,k] - lik.den$R1[l1,k]

      # picking the firms
      lpr_j1    = sj1[I] - lik.den$Rj[l1]
      # if same cluster, needs to remove the PR of the firm in period 1
      if (l1==l2) {
        #lpr_j2    = sj2[I] - log( exp(lik.den$Rj[l2]) - sj1[I])
        lpr_j2     = sj2[I] - lik.den$Rj[l2] - log( 1 - exp(sj1[I]-lik.den$Rj[l2]))
      } else {
        lpr_j2    = sj2[I] - lik.den$Rj[l2]
      }

      # wages
      lpr_y1    = lognormpdf(Y1m[I] , model$A1[l1,k], model$S1[l1,k])
      lpr_y2    = lognormpdf(Y2m[I] , model$A2[l2,k], model$S2[l2,k])

      lpm[I,k] = log(model$pk[k]) + lpr_k1 + lpr_move + lpr_k2 + lpr_y1 + lpr_j1 + lpr_j2 + lpr_y2
    }
  }

  sum(logRowSumExp(lpm))
}

#' computes the total liklihood
#'
#' @param sim jdata and sdata
#' @param model is a model.proba list with A1,S1,A2,S2,pm,g0,g1,pk,nk,nf,Jall_l
#' @export
m2.prob.lik.all <- function(sim,model,ctrl) {
  # ------- UPDATING SOME FIRMS ------- #
  sdata =sim$sdata
  jdata = sim$jdata

  dprior = ctrl$dprior
  model0 = ctrl$model0
  taum   = ctrl$tau

  ### ----- GET MODEL  ---
  nk  = model$nk
  nf  = model$nf
  A1  = model$A1
  S1  = model$S1
  A2  = model$A1
  S2  = model$S1
  g0  = model$g0
  g1  = model$g1
  pk  = model$pk
  pl  = model$pl
  pm  = model$pm

  # ----- GET DATA
  # movers
  Y1m = jdata$y1
  Y2m = jdata$y2
  J1m = jdata$j1
  J2m = jdata$j2
  F1m = jdata$f1
  F2m = jdata$f2
  # stayers
  Y1s = sdata$y1
  Y2s = sdata$y2
  J1s = sdata$j1
  F1s = sdata$f1

  # create the depend variables
  nfirms = sdata[,length(unique(f1))]
  ns = sdata[,.N]
  nm = jdata[,.N]

  # current draw for firm types
  ddtmp  = unique(sdata[,.N,list(f1,j1)])
  Jall_l = ddtmp$j1
  Jall_f = ddtmp$f1
  Jall_s = ddtmp$N

  # extract realized size in cross-section
  sdata[,n1:=.N,f1]
  sizedata = unique(sdata[,list(n1,f1)])
  setkey(sizedata,f1)
  setkey(jdata,f1)
  jdata[,n1:=sizedata[jdata,n1]]
  setkey(jdata,f2)
  jdata[,n2:=sizedata[jdata,n1]]

  N1s = sdata$n1
  N1m = jdata$n1
  N2m = jdata$n2

  ## --- update firms --- #

  ## we prepare a set of constants that will be useful
  ## these are the liklihood in the current (k_1....k_J) for all firms
  dd     = list(Y1s=Y1s,J1s=J1s,J1m=J1m,J2m=J2m,Y1m=Y1m,Y2m=Y2m,F1s=F1s,F1m=F1m,F2m=F2m,N1m=N1m,N2m=N2m,N1s=N1s)
  likden = m2.proba.likden(model,Jall_s,Jall_l)
  rr = m2.prob.lik(dd,model,likden,dprior=dprior)
  return(rr)
}

#' computes the liklihood of movers and stayers
#'
#' @param dd is a list with Y1s, J1s, Y1m, Y2m, J1m, J2m, Jall_f
#' @param model is list with A1,S1,A2,S2,pm,g0,g1,pk,nk,nf,Jall_l
#' @export
m2.prob.lik <- function(dd,model,likden,dprior=1.01) {

  nk = model$nk
  nf = model$nf
  ns = length(dd$Y1s)
  nm = length(dd$Y1m)

  lik_stayers_l = array(0,nf)
  lik_movers_ll = array(0,c(nf,nf))

  ## --- STAYERS --- #
  lps = array(0,c(ns,nk))
  for (l1 in 1:nf) {
    I = which( (dd$J1s==l1))
    if (length(I)==0) next;
    for (k in 1:nk) {
      lpr_y1 = lognormpdf(dd$Y1s[I] , model$A1[l1,k], model$S1[l1,k])  # the wage
      lpr_k1 = model$g0[l1,k] - likden$R0[k]
      lpr_j1 = log(dd$N1s[I]) - likden$Rj[l1]
      lps[I,k] = log(model$pk[k]) + lpr_k1+ lpr_j1 + lpr_y1 + log(1-model$pm[l1,k])
    }
    lik_stayers_l[l1] = sum(logRowSumExp(lps[I,]))
  }
  liks     = sum(logRowSumExp(lps))
  taus     = exp(lps - spread(logRowSumExp(lps),2,nk)) # normalize the k probabilities Pr(k|Y1,Y2,Y3,Y4,l)

  # prior
  liks_prior = (dprior-1) * sum(log(model$pk))
  liks = liks

  ## --- MOVERS --- #
  lpm = array(0,c(nm,nk))
  for (l1 in 1:nf) for (l2 in 1:nf) {
    I = which( (dd$J1m==l1) & (dd$J2m==l2))
    if (length(I)==0) next;
    for (k in 1:nk) {
      lpr_y1    = lognormpdf(dd$Y1m[I] , model$A1[l1,k], model$S1[l1,k])
      lpr_k1    = model$g0[l1,k] - likden$R0[k]
      lpr_j1    = log(dd$N1m[I]) - likden$Rj[l1]
      lpr_move  = log(model$pm[l1,k])
      lpr_k2    = model$g1[l1,l2,k] - likden$R1[l1,k]
      if (l1==l2) {
        lpr_j2    = log(dd$N2m[I]) - likden$Rj[l2] - log( 1 - exp(log(dd$N1m[I])-likden$Rj[l2]))
        if (any(!is.finite(lpr_j2))) browser();
      }  else {
        lpr_j2    = log(dd$N2m[I]) - likden$Rj[l2]
      }
      lpr_y2    = lognormpdf(dd$Y2m[I] , model$A2[l2,k], model$S2[l2,k])

      # sum the log of all
      lpm[I,k] = log(model$pk[k]) + lpr_y1 + lpr_j1 +
        lpr_k1 + lpr_move + lpr_k2 + lpr_j2 + lpr_y2
    }
    lik_movers_ll[l1,l2] = sum(logRowSumExp(lpm[I,]))
  }
  likm     = sum(logRowSumExp(lpm))
  taum     = exp(lpm - spread(logRowSumExp(lpm),2,nk)) # normalize the k probabilities Pr(k|Y1,Y2,Y3,Y4,l)

  # prior
  likm_prior = (dprior-1) * sum(log(model$pk))
  likm = likm

  lik = likm + liks + likm_prior + liks_prior
  return(list(lik=lik,liks=liks,likm=likm,taum=taum,taus=taus,lik_stayers_l=lik_stayers_l,lik_movers_ll=lik_movers_ll))
}

#' precompute the likelihood denominator constants
#' @export
m2.proba.likden <- function(model,Jall_s,Jall_l) {

  # @fixme strcitly speaking the constants in period 1 should also include the movers since it is unconditional on moving.

  ## prepapre denominator constants for cluster choices
  R0 = array(0,c(model$nk))
  for (k in 1:model$nk) {
    R0[k] = logsumexp(model$g0[,k])
  }
  R1 = array(0,c(model$nf,model$nk))
  for (k in 1:model$nk) for (l1 in 1:model$nf) {
    R1[l1,k] = logsumexp(model$g1[l1,,k])
  }
  ## prepapre denominator constants for firm choice within cluster
  Rj = array(0,c(model$nf))
  for (l1 in 1:model$nf) {
    Rj[l1] = logsumexp(log(Jall_s[Jall_l==l1]))
  }
  ## prepapre denominator constants for firm choice within cluster for movers
  #Rjm = array(0,c(model$nf))
  #for (l2 in 1:model$nf) {
  #  Rjm[l2] = logsumexp(log(Jall_s2[Jall_l2==l2]))
  #}

  #return(list(R0=R0,R1=R1,Rj=Rj,Rjm=Rjm))
  return(list(R0=R0,R1=R1,Rj=Rj))
}

#' update the classification using probabilistic model
#'
#' We take a current classification, and value of the model. We then
#' draw new firm classes according to a full probabilistic model.
#' @export
m2.proba.updateclasses <- function(sim,model0,pl,ctrl,debug=FALSE,mobility_smooth=0.0001) {

  # what we need here is to be able to update the classification according
  # to the probabilistic model. This means that I need to compute the posterior
  # probabilities for each firm and each type.

  # this then firest requires to map the BLM model into the probabilistic model.
  model    = m2.proba.new.from.mixt(model0,mobility_smooth=mobility_smooth)
  model$pl = pl
  # we aslo need to extract firm size info (required in the probabilistic model)

  # ------- UPDATING SOME FIRMS ------- #
  sdata =sim$sdata
  jdata = sim$jdata

  dprior = ctrl$dprior
  model0 = ctrl$model0
  taum   = ctrl$tau

  ### ----- GET MODEL  ---
  nk  = model$nk
  nf  = model$nf
  A1  = model$A1
  S1  = model$S1
  A2  = model$A1
  S2  = model$S1
  g0  = model$g0
  g1  = model$g1
  pk  = model$pk
  pl  = model$pl
  pm  = model$pm

  # ----- GET DATA
  # movers
  Y1m = jdata$y1
  Y2m = jdata$y2
  J1m = jdata$j1
  J2m = jdata$j2
  F1m = jdata$f1
  F2m = jdata$f2
  # stayers
  Y1s = sdata$y1
  Y2s = sdata$y2
  J1s = sdata$j1
  F1s = sdata$f1

  # create the depend variables
  nfirms = sdata[,length(unique(f1))]
  ns = sdata[,.N]
  nm = jdata[,.N]

  # current draw for firm types
  ddtmp  = unique(sdata[,.N,list(f1,j1)])
  Jall_l = ddtmp$j1
  Jall_f = ddtmp$f1
  Jall_s = ddtmp$N

  # extract realized size in cross-section
  sdata[,n1:=.N,f1]
  sizedata = unique(sdata[,list(n1,f1)])
  setkey(sizedata,f1)
  setkey(jdata,f1)
  jdata[,n1:=sizedata[jdata,n1]]
  setkey(jdata,f2)
  jdata[,n2:=sizedata[jdata,n1]]

  N1s = sdata$n1
  N1m = jdata$n1
  N2m = jdata$n2
  
  # extract lieklihood of this classification in this model using pl
  lik_before_firm_update_gpr = sim$sdata[,length(unique(f1)),j1][,sum(V1*log(pl[j1]))]

  ## --- update firms --- #

  ## we prepare a set of constants that will be useful
  ## these are the liklihood in the current (k_1....k_J) for all firms
  dd     = list(Y1s=Y1s,J1s=J1s,J1m=J1m,J2m=J2m,Y1m=Y1m,Y2m=Y2m,F1s=F1s,F1m=F1m,F2m=F2m,N1m=N1m,N2m=N2m,N1s=N1s)
  likden = m2.proba.likden(model,Jall_s,Jall_l)
  rr = m2.prob.lik(dd,model,likden,dprior=dprior)
  lik_before_firm_update         = rr$lik
  lik_before_firm_update_movers  = rr$likm
  lik_before_firm_update_stayers = rr$liks
  posterior_pr = array(0,c(ctrl$nfirms_to_update,nf))

  all_res = data.frame()

  # we draw some firms to update
  firms_to_update_indices   = sample.int(nfirms,ctrl$nfirms_to_update)
  firms_to_update_new_types = rep(0,ctrl$nfirms_to_update)
  for (ii in 1:(ctrl$nfirms_to_update)) {
    findex  = firms_to_update_indices[[ii]]
    fidcur  = Jall_f[[findex]]
    lprs    = rep(0,nf) # likelihood stayers
    lprm    = rep(0,nf) # likelihood movers
    lpr     = rep(0,nf) # likelihood total

    l1cur = Jall_l[findex]          # the current class of the firm

    # next, we loop other clusters
    for (l1 in 1:nf) {
      # update constants for the computation of the likelihood
      likden_guess = copy(likden)
      likden_guess$Rj[l1cur]  = log(exp(likden_guess$Rj[l1cur]) - Jall_s[findex])
      likden_guess$Rj[l1]     = log(exp(likden_guess$Rj[l1])    + Jall_s[findex])
      ll_affected = union(l1,l1cur)

      # ------------  STAYERS ------------ #
      # finding all stayers in clusters (l1_cur,l1_guess)
      Is_affected = which( (J1s == l1) | (J1s == l1cur))
      if (ctrl$proba_include[2]==0) Is_affected = which(F1s == fidcur);
      J1s_guess = J1s[Is_affected]

      # update the class for the current firm
      J1s_guess[ F1s[Is_affected] == fidcur ] = l1 # update workers in cur frim

      # compute the likelihood for this group
      lprs[l1] = (ctrl$proba_include[1]==TRUE)*m2.proba.lik.stayers(Y1s[Is_affected],J1s_guess,log(N1s[Is_affected]),model,likden_guess)
      # add the likelihood from the other stayers
      lprs[l1] = lprs[l1] + (ctrl$proba_include[2]==TRUE)*(rr$liks - sum(rr$lik_stayers[ll_affected]))

      # ------------  MOVERS ------------ #
      # finding all movers coming and going to clsuters (l1_cur,l1_guess) and
      Im_affected = which( (J1m == l1) | (J1m == l1cur) | (J2m == l1) | (J2m == l1cur))
      if (ctrl$proba_include[4]==0) Im_affected = which( (F1m == fidcur) | ( F2m== fidcur));
      J1m_guess = J1m[Im_affected]
      J2m_guess = J2m[Im_affected]

      # update the class for the current firm
      J1m_guess[ F1m[Im_affected] == fidcur ] = l1 # update workers moving from current firm
      J2m_guess[ F2m[Im_affected] == fidcur ] = l1 # update workers moving to current firm

      # compute the likelihood
      lprm[l1] = lprm[l1] + (ctrl$proba_include[3]==TRUE)* m2.proba.lik.movers(Y1m[Im_affected],Y2m[Im_affected],J1m_guess,J2m_guess,log(N1m[Im_affected]),log(N2m[Im_affected]),model,likden_guess)
      lprm[l1] = lprm[l1] + (ctrl$proba_include[4]==TRUE)*(rr$likm - sum(rr$lik_movers_ll[,ll_affected]) - sum(rr$lik_movers_ll[ll_affected,]) +
        sum(rr$lik_movers_ll[ll_affected,ll_affected])) # the last is to not remove the intersection twice

      # combine the two
      lpr[l1] = lprs[l1] + lprm[l1]
      posterior_pr[ii,l1] = lpr[l1]
    }

    lpr = lpr + log(pl) # adding the probability of each class among firms
    rownames(posterior_pr) = Jall_f[firms_to_update_indices]
    firms_to_update_new_types[[ii]] = sample.int(nf,1,prob=exp(lpr-max(lpr)))
    if (debug==TRUE) {
      class_true  = sdata[f1==Jall_f[[findex]]][,j1true[1]]
      class_map   = which.max(exp(lpr-max(lpr)))
      class_pre   = Jall_l[[findex]]
      pr_map      = max( exp(lpr-max(lpr))/(sum(exp(lpr-max(lpr))))   )
      str1 = ifelse(class_true==firms_to_update_new_types[[ii]],"*","")
      str2 = ifelse(class_pre!=firms_to_update_new_types[[ii]],"!","")
      flog.info("j1true=%2i j1pre=%2i j1map=%2i j1new=%2i pr_map=%2f %s%s",class_true,class_pre,class_map, firms_to_update_new_types[[ii]],pr_map,str2,str1)
      tmp=0
      all_res = rbind(all_res,data.frame(f1=fidcur,l1=1:model$nf,liks=lprs,likm=lprm,l1pre=Jall_l[[findex]],l1true=class_true))
    } else {
      class_map   = which.max(exp(lpr-max(lpr)))
      class_pre   = Jall_l[[findex]]
      pr_map      = max( exp(lpr-max(lpr))/(sum(exp(lpr-max(lpr))))   )
      str2 = ifelse(class_pre!=firms_to_update_new_types[[ii]],"!","")
      flog.info("j1pre=%2i j1map=%2i j1new=%2i pr_map=%2f %s",class_pre,class_map, firms_to_update_new_types[[ii]],pr_map,str2)
      all_res = rbind(all_res,data.frame(f1=fidcur,l1=1:model$nf,liks=lprs,likm=lprm,l1pre=Jall_l[[findex]]))
    }

    if (ii%%50==0) flog.info("computing firms posterior %i/%i",ii,ctrl$nfirms_to_update);
  }
  # construct the updated set of firms ID
  Jall_l_new = copy(Jall_l)
  for (ii in 1:(ctrl$nfirms_to_update)) {
    findex = firms_to_update_indices[[ii]] # index of the firm to update
    Jall_l_new[[findex]] = firms_to_update_new_types[[ii]]
  }
  names(Jall_l_new) = Jall_f


  flog.info("done updating firms (%i/%i changed)", sum(Jall_l_new!=Jall_l), ctrl$nfirms_to_update)

  return(list(grps=Jall_l_new,all=all_res,updates=sum(Jall_l_new!=Jall_l),
              lik_before_firm_update=lik_before_firm_update,posterior_pr=posterior_pr,lik_before_firm_update_movers=lik_before_firm_update_movers,lik_before_firm_update_stayers=lik_before_firm_update_stayers,lik_before_firm_update_gpr=lik_before_firm_update_gpr))
}

#' computes simple posterior probabilities
#' @export
m2.prob.get_firm_posteriors <- function(sim,model) {

  nk = model$nk
  nf = model$nf

  lik_stayers_l = array(0,nf)
  lik_movers_ll = array(0,c(nf,nf))

  R0 = array(0,c(model$nk))
  for (k in 1:model$nk) {
    R0[k] = logsumexp(model$g0[,k])
  }

  ## we can just loop over l1 and k
  Y1s = sim$sdata$y1
  ns = sim$sdata[,.N]
  dd = sim$sdata[,list(f1)]
  nfids = length(unique(sim$sdata$f1))
  lpr_post = array(0,c(nfids,nf))

  ## --- compute posterior for each i,l1 --- #
  lps      = array(0,c(ns,nk))
  lps_post = array(0,c(ns,nf))
  for (l1 in 1:nf) {
    for (k in 1:nk) {
      lpr_y1 = lognormpdf(Y1s , model$A1[l1,k], model$S1[l1,k])  # the wage
      lpr_k1 = model$g0[l1,k] - R0[k]
      lps[,k] = log(model$pk[k]) + lpr_k1 + lpr_y1 + log(1-model$pm[l1,k])
    }
    lps_post[,l1]     = logRowSumExp(lps)

    # aggregate across firms
    dd[, lpost := lps_post[,l1]]
    lpr_post[,l1] = dd[,sum(lpost),f1][,V1]
  }

  rownames(lpr_post) = dd[,sum(lpost),f1][,f1]
  #lpr_post = lpr_post - spread(logRowSumExp(lpr_post),2,nf)

  return(lpr_post)
}

#' precompute the likelihood denominator constants
#' @export
m2.proba.importancelik <- function(sim,model,eps_mix=0.1,class_draws=200,mobility_smooth=0.0001) {

  fids = sim$sdata[,unique(f1)]
  ctrl      = em.control(tol=1e-7,check_lik=F,dprior=1.01,fixb=T)
  
  flog.info("computing posterior probabilities for all %i firms",length(fids))
  model.proba   = m2.proba.new.from.mixt(model,mobility_smooth=mobility_smooth)
  posterior_pr = m2.prob.get_firm_posteriors(sim,model.proba)
  class_pr      = model$NNs + rowSums(model$NNm) # @fixme, this should be frequency at the firm level?
  class_pr      = as.numeric(class_pr / sum(class_pr))
  
  # we draw classifications, evaluate the likelihood
  rr = data.frame()
  posterior_pr = posterior_pr - logRowSumExp(posterior_pr)
  lpr          = log( eps_mix * spread(class_pr,1,length(fids)) + (1-eps_mix) *  exp(posterior_pr))

  flog.info("drawing %i classifications (eps=%f)",class_draws,eps_mix)
  for (i in 1:class_draws) {

    grps_class = rep(0,length(fids))
    names(grps_class) = rownames(posterior_pr)

    # draw a classification from the posterior probabilities
    # and compute classification prior probability
    log_proposal_pr = 0
    log_classificiation_prior = 0
    for (j in 1:length(fids)) {
      l1 = sample.int(model$nf,1,prob=exp(lpr[j,]))
      grps_class[j] = l1
      log_proposal_pr = log_proposal_pr + as.numeric(lpr[j,l1])
      log_classificiation_prior = log_classificiation_prior + log(class_pr[l1])
    }

    # assign the classification to the data
    sim            = grouping.append(sim,grps_class)

    # compute the likelihood
    res = m2.prob.lik.all(sim,model.proba,ctrl)

    rr = rbind(rr,data.frame(lik=res$lik,prop_pr = log_proposal_pr,class_prior=log_classificiation_prior,eps=eps_mix))
    if (i%%20==0) flog.info("drawing classifictions %i/%i",i,class_draws)
  }

  rr = data.table(rr[is.finite(rr$lik),])

  res = list()
  res$model        = model
  res$model.proba  = model.proba
  res$posterior_pr = posterior_pr
  res$all          = rr
  res$eps_mix      = eps_mix

  res$lik = rr[is.finite(lik), logsumexp(lik - prop_pr + class_prior)]

  res
}

#' precompute the likelihood denominator constants
#' @export
m2.proba.gibbs <- function(sim,model,pl,eps_mix=0.1,mobility_smooth=1e-10,nfirms_to_update=1,maxiter=100,grp_start=NA,appendto=NA,debug=F) {
  
  fids = sim$sdata[,unique(f1)]
  ctrl      = em.control(tol=1e-7,check_lik=F,dprior=1.01,fixb=T)
  ctrl$proba_include =c(1,1,1,1) #use full likelihood to upate
  ctrl$nfirms_to_update = nfirms_to_update
  
  # USING RANDOM CLASSIFICATION AS STARTING POINT
  if (any(is.na(grp_start))) {
    grps_random = sample(model$nf,length(fids),replace=T)
    names(grps_random) = fids
    sim            = grouping.append(sim,grps_random)
  } else {
    sim            = grouping.append(sim,grp_start)
  }

  # check if we should append
  if (is.na(appendto)) {
    rr = data.frame()
    istart=1
  } else{
    rr = appendto
    istart  = max(appendto$i)+1
    maxiter = istart + maxiter
  }
  
  # run the GIBBS SAMPLER
  for (i in istart:maxiter) {
    grps_r         = m2.proba.updateclasses(sim,model,pl,ctrl,debug=debug,mobility_smooth=mobility_smooth)
    sim            = grouping.append(sim,grps_r$grps,sort=FALSE)
    #print(AA)
    if (debug==T) {
      AA = acast(sim$sdata[,.N,list(j1,j1true)],j1~j1true,fill=0,value.var = "N")
      rr = rbind(rr,data.frame(i=i,step = i * ctrl$nfirms_to_update, match = sum(diag(AA))/sum(AA), mse = sim$sdata[,mean( (j1-j1true)^2)], fsize= fsize, lik_pre=grps_r$lik_before_firm_update))
      flog.info("iter=%i mse=%f match=%f lik_pre=%4.4f",i * ctrl$nfirms_to_update,sim$sdata[,mean( (j1-j1true)^2)],sum(diag(AA))/sum(AA),grps_r$lik_before_firm_update)
    } else {
      rr = rbind(rr,data.frame(i=i,step = i * ctrl$nfirms_to_update, lik=grps_r$lik_before_firm_update, likm=grps_r$lik_before_firm_update_movers,liks=grps_r$lik_before_firm_update_stayers,likg=grps_r$lik_before_firm_update_gpr,updates=grps_r$updates))
      flog.info("iter=%i lik_pre=%4.4e movers=%4.4e stayers=%4.4e grps=%4.4e",i * ctrl$nfirms_to_update,grps_r$lik_before_firm_update,grps_r$lik_before_firm_update_movers,grps_r$lik_before_firm_update_stayers,grps_r$lik_before_firm_update_gpr)
    }
  }
  
  return(list(rr=rr,last_grp = grps_r$grps))
}

#' extract the proportions from a classification
#' @export
m2.proba.getlcasspr <- function(sim,nf) {
  
  V = rep(0,nf)
  # get the number of firms per class 
  for (i in 1:nf) {
    V[i] = sim$sdata[j1==i,length(unique(f1))]
  }
  class_pr = V/sum(V)
  lik = sum(log(class_pr)*V)
  
  return(list(class_pr = class_pr,lik=lik,NNf=V))
}
