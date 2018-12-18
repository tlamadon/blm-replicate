# Creates, simulates and estimates
# the 4 period interactive model

#' Creates a 2 period mini model. Can be set to be linear or not, with serial correlation or not
#' @export
m4.mini.new <- function(nf,linear=FALSE,serial=F,r1=0.3 + 0.6*runif(1) ,r4=0.3 + 0.6*runif(1),fixb=F,eps_sd_factor=1) {
  m = list()
  m$nf = nf

  # generating intercepts and interaction terms 
  m$A1  = rnorm(nf)
  m$B1  = exp(rnorm(nf)/2) 
  m$A2  = rnorm(nf)
  m$A2s = rnorm(nf)
  m$B2  = exp(rnorm(nf)/2) 
  m$A3  = rnorm(nf)
  m$A3s = rnorm(nf)
  m$B3  = exp(rnorm(nf)/2) 
  m$A4  = rnorm(nf)
  m$B4  = exp(rnorm(nf)/2) 

  if (fixb) {
    m$B2=m$B1
    m$B3=m$B1
    m$B4=m$B1
  }
  
  # generat the effect of the future firm
  m$C2 = rnorm(nf)
  m$C3 = rnorm(nf)

  # generat the auto-cor
  m$r1 = r1
  m$r4 = r4

  # generate the E(alpha|l,l') and sd
  m$EEm  = array(rnorm(nf^2),c(nf,nf))
  m$EEsd = exp(array(rnorm(nf^2),c(nf,nf))/2)

  # generate the variance of eps
  m$eps2s_sd  = exp(rnorm(nf)/2) * eps_sd_factor
  m$eps3s_sd  = exp(rnorm(nf)/2) * eps_sd_factor
  #m$eps2m_sd  = array(exp(rnorm(nf*nf)/2),c(nf,nf)) * eps_sd_factor
  #m$eps3m_sd  = array(exp(rnorm(nf*nf)/2),c(nf,nf)) * eps_sd_factor
  m$eps2m_sd  = spread(array(exp(rnorm(nf)/2),c(nf)) * eps_sd_factor,2,nf)
  m$eps3m_sd  = spread(array(exp(rnorm(nf)/2),c(nf)) * eps_sd_factor,1,nf)
  m$nu1_sd  = sqrt(  (1-r1^2) *  m$eps2s_sd^2 )
  m$nu2_sd  = sqrt(  (1-r4^2) *  m$eps3s_sd^2 )
  m$eps_cor = runif(1)

  # generate the E(alpha|l) and sd
  m$Em   = sort(rnorm(nf))
  m$Esd  = exp(rnorm(nf)/2)

  m$rho32m = 0.2*0.5*runif(1)
  m$rho32s = 0.2*0.5*runif(1)
  
  # set normalization
  rr = m$B1[1] - m$r1 * m$B2[1]
  m$B1 = m$B1/rr
  m$B2 = m$B2/rr
  m$B3 = m$B3/rr
  m$B4 = m$B4/rr
  m$A4[nf] = m$r4 * m$A3[nf]
  m$C2[1]  = 0
  m$C3[1]  = 0
  m$A1[1] = m$r1*m$A2[1]
  
  m$Nm = array(5000 ,c(m$nf,m$nf)) + diag(rep(5000,m$nf))
  m$Ns = array(20000,m$nf)
  
  if (linear) {
    m$B1[]  = 1
    m$B2[]  = 1
    m$B3[]  = 1
    m$B4[]  = 1
  }
  
  return(m)
}

#' simulate movers according to the model
#' @export
m4.mini.simulate.all <- function(model,N=50000) {
  
  nf = model$nf
  D   = sample.int(2*nf,N,replace=T,prob = c(model$Ns,rowSums(model$Nm)))
  Mb  = (D>nf)*1
  J1  = D - Mb*nf
  J2  = rep(0,N)
  rr  = data.table(j1=J1,m=Mb,j2=as.integer(0))
  rr[m==1, j2:= sample.int(nf,.N,replace=T,prob=model$Nm[j1,]) ,j1]
  Ns  = rr[m==0,.N,j1][order(j1),N]
  Nm  = acast(rr[m==1,.N,list(j1,j2)],j1~j2,fill=0,value.var = "N")
  
  jdata = m4.mini.simulate.movers(model,Nm)
  sdata = m4.mini.simulate.stayers(model,Ns)
  jdata[,m:=1]
  sdata[,m:=0]
  return(rbind(jdata[,list(alpha,y1,y2,y3,y4,j1,j2,m)],sdata[,list(alpha,y1,y2,y3,y4,j1,j2=j1,m)])) 
}



#' simulate movers according to the model
#' @export
m4.mini.simulate.movers <- function(model,NNm) {

  J1 = array(0,sum(NNm)) 
  J2 = array(0,sum(NNm)) 
  Y1 = array(0,sum(NNm)) 
  Y2 = array(0,sum(NNm)) 
  Y3 = array(0,sum(NNm)) 
  Y4 = array(0,sum(NNm)) 
  e2 = array(0,sum(NNm)) 
  e3 = array(0,sum(NNm)) 
  K  = array(0,sum(NNm)) 
  
  A1  = model$A1
  B1  = model$B1
  A2  = model$A2
  B2  = model$B2
  A3  = model$A3
  B3  = model$B3
  A4  = model$A4
  B4  = model$B4
  C2  = model$C2
  C3  = model$C3
  EEm = model$EEm
  EEsd= model$EEsd
  eps2m_sd=model$eps2m_sd
  eps3m_sd=model$eps3m_sd
  nu1_sd=model$nu1_sd
  nu2_sd=model$nu2_sd
  rho32m = model$rho32m
  nf  = model$nf
  r1 = model$r1
  r4 = model$r4
  
  i =1 
  for (l1 in 1:nf) for (l2 in 1:nf) {
    if (NNm[l1,l2]==0) next;
    I = i:(i+NNm[l1,l2]-1)
    ni = length(I)
    jj = l1 + nf*(l2 -1)
    J1[I] = l1
    J2[I] = l2
    
    # draw alpha
    K[I] = EEm[l1,l2] + EEsd[l1,l2]*rnorm(ni)
    
    # draw epsilon2 and epsilon3 corolated
    eps1 = eps2m_sd[l1,l2] * rnorm(ni)
    if (eps2m_sd[l1,l2]>0) {
      eps2 = eps3m_sd[l1,l2]/eps2m_sd[l1,l2]*rho32m*eps1    +    sqrt(  (1-rho32m^2) *  eps3m_sd[l1,l2]^2 ) * rnorm(ni)
    } else {
      eps2 = sqrt(  (1-rho32m^2) *  eps3m_sd[l1,l2]^2 ) * rnorm(ni)
    }
    # compute Y2 and Y3
    Y2[I]  = A2[l1] + B2[l1]*K[I] + C2[l2] + eps1 
    Y3[I]  = A3[l2] + B3[l2]*K[I] + C3[l1] + eps2 
    e2[I]  = eps1
    e3[I]  = eps2

    # compute Y1,Y4
    Y1[I]  = r1*Y2[I] + A1[l1] - r1*A2[l1] + (B1[l1] - r1 * B2[l1])*K[I] +  nu1_sd[l1] * rnorm(ni)
    Y4[I]  = r4*Y3[I] + A4[l2] - r4*A3[l2] + (B4[l2] - r4 * B3[l2])*K[I] +  nu2_sd[l2] * rnorm(ni)
        
    i = i + NNm[l1,l2]
  }
  
  jdatae = data.table(alpha=K,y1=Y1,y2=Y2,y3=Y3,y4=Y4,j1=J1,j2=J2,e2=e2,e3=e3)
  return(jdatae)   
}

#' simulate statyers according to the model
#' @export
m4.mini.simulate.stayers <- function(model,NNs,ks=c(1,1)) {
  
  J1 = array(0,sum(NNs)) 
  Y1 = array(0,sum(NNs)) 
  Y2 = array(0,sum(NNs)) 
  Y3 = array(0,sum(NNs)) 
  Y4 = array(0,sum(NNs)) 
  Eps2 = array(0,sum(NNs)) 
  Eps3 = array(0,sum(NNs)) 
  Nu1 = array(0,sum(NNs)) 
  Nu4 = array(0,sum(NNs)) 

  K  = array(0,sum(NNs)) 
  
  A1   = model$A1
  B1   = model$B1
  A2s  = model$A2s
  A2   = model$A2
  B2   = model$B2
  A3   = model$A3
  A3s  = model$A3s
  B3   = model$B3
  A4   = model$A4
  B4   = model$B4
  r1   = model$r1
  r4   = model$r4
  Em   = model$Em
  Esd  = model$Esd
  eps2s_sd=model$eps2s_sd
  eps3s_sd=model$eps3s_sd
  rho32s=model$rho32s
  nu1_sd=model$nu1_sd
  nu2_sd=model$nu2_sd
  
  nf  = model$nf
  
  i =1 
  for (l1 in 1:nf) {
    I = i:(i+NNs[l1]-1)
    ni = length(I)
    J1[I] = l1
    
    # draw alpha
    K[I] = Em[l1] + Esd[l1]*rnorm(ni)
    
    # draw epsilon2 and epsilon3 corolated
    eps1 = eps2s_sd[l1] * rnorms(ni,ks)
    if (eps2s_sd[l1]>0) {
      eps2 = eps3s_sd[l1]/eps2s_sd[l1]*rho32s*eps1    +    sqrt(  (1-rho32s^2) *  eps3s_sd[l1]^2 ) * rnorms(ni,ks)
    } else {
      eps2 = sqrt(  (1-rho32s^2) *  eps3s_sd[l1]^2 ) * rnorms(ni,ks)
    }
    # compute Y2 and Y3
    Y2[I]  = A2s[l1] + B2[l1]*K[I] + eps1 
    Y3[I]  = A3s[l1] + B3[l1]*K[I] + eps2 
    
    # compute Y1,Y4
    Nu1[I] = nu1_sd[l1] * rnorms(ni,ks)
    Nu4[I] = nu2_sd[l1] * rnorms(ni,ks)
    Y1[I]  = r1*Y2[I] + A1[l1] - r1*A2[l1] + (B1[l1] - r1 * B2[l1])*K[I] + Nu1[I] # here same as movers
    Y4[I]  = r4*Y3[I] + A4[l1] - r4*A3[l1] + (B4[l1] - r4 * B3[l1])*K[I] + Nu4[I] # here same as movers
    #Y1[I]  = A1[l1] + B1[l1]*K[I]     + r1 * (Y2[I] - A2[l1] - B2[l1]*K[I])    + Nu1[I]
    #Y4[I]  = r4*Y3[I] + A4[l1] - r4*A3[l1] + (B4[l1] - r4 * B3[l1])*K[I] + Nu4[I]
    
    Eps2[I]=eps1
    Eps3[I]=eps2

    i = i + NNs[l1]
  }
  
  sdatae = data.table(alpha=K,y1=Y1,y2=Y2,y3=Y3,y4=Y4,j1=J1,eps2=Eps2,eps3=Eps3,nu1=Nu1,nu4=Nu4)
  return(sdatae)  
}

#' compute the posterior probabilities of staying and moving to any
#' of the j2 firms for a given (j1,a,y2)
#' @export
m4.mini.posteriormobility <- function(model,j1,a,y2,include_y2=1,include_a=1) {
  # Pr[Y2_i|m,j1_i,j2,alpha_i] for m=0 and j2=1..10
  #p1 = c( (Y2[i] - model$A2s[l1] - model$B2[l1]*a)/model$eps2s_sd[l1], 
  #        (Y2[i] - model$A2[l1]  - model$B2[l1]*a - model$C2)/ model$eps2m_sd[l1,])     
  p1 = lognormpdf(y2, 
                  c(model$A2s[j1] + model$B2[j1]*a,model$A2[j1] + model$B2[j1]*a + model$C2),
                  c(model$eps2s_sd[j1],model$eps2m_sd[j1,]))
  
  # Pr[alpha_i|m,j1_i,j2]
  #p2 = c( (a - model$Em[j1])/model$Esd[j1],
  #        (a - model$EEm[j1,])/model$EEsd[j1,])
  #p2 = lognormpdf(p2)
  p2 = lognormpdf(a,c(model$Em[j1],model$EEm[j1,]), c(model$Esd[j1],model$EEsd[j1,]))
  
  # Pr[m=1,j2| j1_i]
  p3 = c(model$Ns[j1],model$Nm[j1,])
  p3 = p3/sum(p3)
  p3 = log(p3)
  
  # combine the probabilities
  pp = p3 + p2*include_a + p1*include_y2
  #pp = p3 + p2 + p1
  pp = pp - logsumexp(pp)
  
  return(pp)
}

#' compute E(y2 | alpha, j1), alpha should be a vector
#' @export
m4.mini.y2_cond_aj1 <- function(model,alpha,j1) {
  
  N   = length(alpha)
  nf  = model$nf
  
  # get the PR[m,j2|j1]
  pm = c(model$Ns[j1],model$Nm[j1,])
  lpm = spread(log(as.numeric(pm/sum(pm))),1,N)
  
  # get the Pr[alpha|m,j1,j2]
  lpa = lognormpdf( spread(alpha,2,nf+1) ,
                    spread( c(model$Em[j1], model$EEm[j1,] ), 1, N), 
                    spread( c(model$Esd[j1],model$EEsd[j1,]), 1, N))
  
  # get the Pr[j2,m|alpha,j1]
  lp  = lpa+lpm - spread(blm:::logRowSumExp(lpa+lpm),2,nf+1)
  
  # finally we multiply by E(y2|alpha,j1,j2,m)
  Y2c = c( model$A2s[j1] + model$B1[j1]* alpha, model$A2[j1] + model$B1[j1]*spread(alpha,2,nf) + spread(model$C2,1,N )  )
  Y2 = rowSums( Y2c*exp(lp))/rowSums(exp(lp))
  
  return(Y2)
}

#' compute E(y2 | alpha, j1), alpha should be a vector
#' @export
m4.mini.y2_cond_a <- function(model,alpha) {
  
  N   = length(alpha)
  nf  = model$nf
  
  pr_j1 = c(model$Ns,rowSums(model$Nm))
  pr_j1 = pr_j1/sum(pr_j1)
  Y2 = 0
  
  for (j1 in 1:nf) {
    Y2 = Y2 + pr_j1[j1]*m4.mini.y2_cond_aj1(model,alpha,j1)
  }

  return(Y2)
}



#' generate the wages
#' @export
m4.mini.simulatewages <- function(m,j1,a,y2,move,j2) {
  # stayers
  if (move==0) {
    e2  = y2 - (m$A2s[j1]  + m$B2[j1]*a)
    e3  = m$eps3s_sd[j1]/m$eps2s_sd[j1]*m$rho32s * e2 + sqrt(  (1-m$rho32s^2) *  m$eps3s_sd[j1]^2 ) * rnorm(1)
    y3  =              m$A3s[j1]  + m$B3[j1]*a + e3
    y1  = m$r1*y2    + m$A1[j1] - m$r1*m$A2[j1] + (m$B1[j1] - m$r1 * m$B2[j1])*a + m$nu1_sd[j1] * rnorm(1)
    y4  = m$r4*y3    + m$A4[j1] - m$r4*m$A3[j1] + (m$B4[j1] - m$r4 * m$B3[j1])*a + m$nu2_sd[j1] * rnorm(1)
  } else {
    e2  = y2 - (m$A2[j1] + m$B2[j1]* a + m$C2[j2])
    e3  = m$eps3m_sd[j1,j2]/m$eps2m_sd[j1,j2]*m$rho32m*e2 + sqrt(  (1-m$rho32m^2) *  m$eps3m_sd[j1,j2]^2 ) * rnorm(1)
    y3  =            m$A3[j2] + m$B3[j2]* a + m$C3[j1] + e3
    y1  = m$r1*y2  + m$A1[j1] - m$r1*m$A2[j1] + (m$B1[j1] - m$r1 * m$B2[j1])*a + m$nu1_sd[j1] * rnorm(1)
    y4  = m$r4*y3  + m$A4[j2] - m$r4*m$A3[j2] + (m$B4[j2] - m$r4 * m$B3[j2])*a + m$nu2_sd[j2] * rnorm(1)
  }
  
  return(list(y1=y1,y3=y3,y4=y4))
}

#' Simulates counterfactual
#' @export
m4.mini.simulate.counter <- function(model,N=10000,d_e2=0,include_y2=1,include_a=1,seed=0,d_j1=0,d_alpha=0) {

  model$eps2m_sd  = pmax(model$eps2m_sd,0.001)
  model$EEsd  = pmax(model$EEsd,0.001)
    
  if (seed>0) set.seed(seed);
  
  # -------------- COMPUTEING BASELINE DISTRIBUTION F(alpha,j1,e2) ------------- #
  nf = model$nf
  
  # draw (j1,m) from the joint from the data
  D   = sample.int(2*nf,N,replace=T,prob = c(model$Ns,rowSums(model$Nm)))
  Mb  = (D>nf)*1
  J1b = D - Mb*nf
  J2  = rep(0,N)
  
  # prepare variables
  A     = rep(0,N)  # worker type
  Ab    = rep(0,N)  # worker type in the baseline
  E2    = rep(0,N)  # residual in period 2
  E2b   = rep(0,N)  # residual in period 2 in the baseline
  Y2    = rep(0,N)  # wage in period 2
  J2b   = rep(0,N)  # j2 in the baseline
  P     = rep(0,N) 
  
  # for the stayers
  I      = which(Mb==Mb)
  Ab[I]  = model$Em[J1b[I]] + rnorm(length(I))* model$Esd[J1b[I]]
  E2b[I]  = rnorm(length(I))*model$eps2s_sd[J1b[I]]
  #Ab   = model$Em[J1b] + rnorm(N)* model$Esd[J1b]
  #E2b  = rnorm(N)*model$eps2s_sd[J1b]
  
  # # for the movers, we need to integrate the J2
  for (j1 in 1:nf) {
    I = which((J1b==j1) & (Mb==1))
    if (length(I)==0) next;
    j2    = sample.int(nf,length(I),replace=T,prob = model$Nm[j1,])
    J2b[I] = j2
    Ab[I]  = model$EEm[j1,j2] + rnorm(length(I))*model$EEsd[j1,j2]
    E2b[I] = rnorm(length(I))*model$eps2m_sd[j1,j2]
  }

  # -------------- SETTING THE COUNTERFACTUAL ------------- #

  A  = Ab  + d_alpha               # shift alpha 
  E2 = E2b + d_e2                  # shift epsilon2
  J1 = pmax(pmin(J1b+d_j1,nf),1)    # shift j1
  J2 = J2b
  
  # -------------- SIMULATING THE COUNTERFACTUAL ------------- #
  
  #  ----  Compute Y2 ---- #
  I     = which(Mb==Mb)
  Y2[I] = model$A2s[J1[I]] + model$B2[J1[I]]* A[I] + E2[I]
  for (j1 in 1:nf) {
    I = which((J1==j1) & (Mb==1))
    if (length(I)==0) next;
    Y2[I] = model$A2[j1] + model$B2[j1]*A[I] + model$C2[J2[I]] + E2[I]
  }
  # Y2 = model$A2s[J1] + model$B2[J1]* A + E2
  
  # ----- Compute (m,j2,Y1,Y2,Y4) ---- #
  Y1 = rep(0,N)
  Y3 = rep(0,N)
  Y4 = rep(0,N)
  M  = rep(0,N)
  for (i in 1:N) {
    
    # draw mobility
    pp    = m4.mini.posteriormobility(model,J1[i],A[i],Y2[i],include_y2=include_y2,include_a=include_a)
    P[i]  = 1-exp(pp[1])/sum(exp(pp))
    M[i]  = sample.int(nf+1,1,prob=exp(pp))-1
    J2[i] = M[i] * (M[i]>=1)
    M[i]  = 1 * (M[i]>=1)
    
    # draw wages
    ww = m4.mini.simulatewages(model,J1[i],A[i],Y2[i],M[i],J2[i])
    
    Y1[i] = ww$y1
    Y3[i] = ww$y3
    Y4[i] = ww$y4
  }
  
  dd = data.table(y1=Y1,y2=Y2,y3=Y3,y4=Y4,j1=J1,j2=J2,a=A,m=M,mb=Mb,j2b = J2b,pm=P,e2=E2) 
  return(dd)
}


#' impute movers according to the model
#' @export
m4.mini.impute.movers <- function(model,jdatae) {
    
  A1  = model$A1
  B1  = model$B1
  A2  = model$A2
  B2  = model$B2
  A3  = model$A3
  B3  = model$B3
  A4  = model$A4
  B4  = model$B4
  C2  = model$C2
  C3  = model$C3
  EEm = model$EEm
  EEsd= model$EEsd
  eps2m_sd=model$eps2m_sd
  eps3m_sd=model$eps3m_sd
  nu1_sd=model$nu1_sd
  nu2_sd=model$nu2_sd
  rho32m = model$rho32m
  nf  = model$nf
  r1 = model$r1
  r4 = model$r4
  
  # ------  impute K, Y1, Y4 on jdata ------- #
  jdatae.sim = copy(jdatae)
  jdatae.sim[, c('k_imp','y1_imp','y2_imp','y3_imp','y4_imp','e2','e3') := {
    ni = .N
    jj = j1 + nf*(j2-1)
    Ki = EEm[j1,j2] + EEsd[j1,j2]*rnorm(ni)
 
    # draw epsilon2 and epsilon3 corolated
    eps1 = eps2m_sd[j1,j2] * rnorm(ni)
    if (eps2m_sd[j1,j2]>0) {
      eps2 = eps3m_sd[j1,j2]/eps2m_sd[j1,j2]*rho32m*eps1    +    sqrt(  (1-rho32m^2) *  eps3m_sd[j1,j2]^2 ) * rnorm(ni)
    } else {
      eps2 = sqrt(  (1-rho32m^2) *  eps3m_sd[j1,j2]^2 ) * rnorm(ni)
    }
    Y2  = A2[j1] + B2[j1]*Ki + C2[j2] + eps1 
    Y3  = A3[j2] + B3[j2]*Ki + C3[j1] + eps2 

    # compute Y1,Y4
    Y1  = r1*Y2 + A1[j1] - r1*A2[j1] + (B1[j1] - r1 * B2[j1])*Ki +  nu1_sd[j1] * rnorm(ni)
    Y4  = r4*Y3 + A4[j2] - r4*A3[j2] + (B4[j2] - r4 * B3[j2])*Ki +  nu2_sd[j2] * rnorm(ni)

    list(Ki,Y1,Y2,Y3,Y4,eps1,eps2)
  },list(j1,j2)]

  return(jdatae.sim)   
}

#' impute stayers according to the model
#' @export
m4.mini.impute.stayers <- function(model,sdatae,ks=c(1,1)) {
    
  A1  = model$A1
  B1  = model$B1
  A2s  = model$A2s
  A2  = model$A2
  B2  = model$B2
  A3s  = model$A3s
  A3  = model$A3
  B3  = model$B3
  A4  = model$A4
  B4  = model$B4
  C2  = model$C2
  C3  = model$C3
  Em = model$Em
  Esd= model$Esd
  eps2s_sd=model$eps2s_sd
  eps3s_sd=model$eps3s_sd
  nu1_sd=model$nu1_sd
  nu2_sd=model$nu2_sd
  rho32s = model$rho32s
  nf  = model$nf
  r1 = model$r1
  r4 = model$r4
  
  # ------  impute K, Y1, Y4 on jdata ------- #
  sdatae.sim = data.table(copy(sdatae))
  sdatae.sim[, c('k_imp','y1_imp','y2_imp','y3_imp','y4_imp','e2','e3') := {
    ni = .N
    Ki = Em[j1] + Esd[j1]*rnorm(ni)
 
    # draw epsilon2 and epsilon3 corolated
    eps1 = eps2s_sd[j1] * rnorms(ni,ks)
    if (eps2s_sd[j1]>0) {
      eps2 = eps3s_sd[j1]/eps2s_sd[j1]*rho32s*eps1    +    sqrt(  (1-rho32s^2) *  eps3s_sd[j1]^2 ) * rnorms(ni,ks)
    } else {
      eps2 = sqrt(  (1-rho32s^2) *  eps3s_sd[j1]^2 ) * rnorms(ni,ks)
    }
    Y2  = A2s[j1] + B2[j1]*Ki + eps1 
    Y3  = A3s[j1] + B3[j1]*Ki + eps2 

    # compute Y1,Y4
    Y1  = r1*Y2 + A1[j1] - r1*A2[j1] + (B1[j1] - r1 * B2[j1])*Ki +  sqrt(  (1-r1^2) *  nu1_sd[j1]^2 ) * rnorms(ni,ks)
    Y4  = r4*Y3 + A4[j1] - r4*A3[j1] + (B4[j1] - r4 * B3[j1])*Ki +  sqrt(  (1-r4^2) *  nu2_sd[j1]^2 ) * rnorms(ni,ks)

    list(Ki,Y1,Y2,Y3,Y4,eps1,eps2)
  },list(j1)]

  return(sdatae.sim)   
}

# simulate movers according to the model
test.simulate <- function() {
  jdata[, mean( model$A1[j1] + model$EEm[j1,j2] + model$r1*model$C2[j2] - y1   ),list(j1,j2)]
  jdata[, mean( model$A2[j1] + model$EEm[j1,j2] + model$C2[j2] - y2   ),list(j1,j2)]
  jdata[, mean( model$A3[j2] + model$EEm[j1,j2] + model$C3[j1] - y3   ),list(j1,j2)]
  jdata[, mean( model$A4[j2] + model$EEm[j1,j2] + model$r4*model$C3[j1] - y4   ),list(j1,j2)]
  
  model = m4.mini.new(10,linear = T)
  NNs  = array(100000/model$nf,model$nf)
  sdata = m4.mini.simulate.stayers(model,NNs)
  
  # checking stayers
  sdata[, mean(  alpha - model$Em[j1] ),list(j1)]
  sdata[, mean( model$A2s[j1] + model$Em[j1]  - y2   ),list(j1)]
  sdata[, with(model,mean( r1*y2 + A1[j1]-r1*A2[j1] + (B1[j1]-r1*B2[j1])*Em[j1]  - y1   )),list(j1)]
  sdata[, with(model,mean( r4*y3 + A4[j1]-r4*A3[j1] + (B4[j1]-r4*B3[j1])*Em[j1]  - y4   )),list(j1)]
  
  jdata[, mean( model$A1[j1] + model$EEm[j1,j2] + model$r1*model$C2[j2] - y1   ),list(j1,j2)]
  jdata[, mean( model$A2[j1] + model$EEm[j1,j2] + model$C2[j2] - y2   ),list(j1,j2)]
  jdata[, mean( model$A3[j2] + model$EEm[j1,j2] + model$C3[j1] - y3   ),list(j1,j2)]
  jdata[, mean( model$A4[j2] + model$EEm[j1,j2] + model$r4*model$C3[j1] - y4   ),list(j1,j2)]
}


# ==================================  ESTIMATION =============================

#' this function uses LIML to estimate a model with a given normalization
#' and we want to do it in a non-stationary way
m4.mini.liml.int <- function(Y1,Y2,J1,J2,norm=1,fixb=F,r1=0,r4=0,alpha=0) {
    
  L = max(J1)
  N = length(Y1)
  
  # ----- prepare the different regressors ----- #
  D1 = array(0,c(N,L))
  I = (1:N) + N*(J1-1)
  D1[I]=1
  D2 = array(0,c(N,L))
  I = (1:N) + N*(J2-1)
  D2[I]=1

  # contrust full matrix
  if (fixb) {
    X1 = D1*spread(Y1,2,L)/(1-r1)-D2*spread(Y2,2,L)/(1-r4)  # construct the matrix for the interaction terms    
    xdim = L
  } else {
    X1 = cbind(D1*spread(Y1,2,L),D2*spread(Y2,2,L))  # construct the matrix for the interaction terms    
    xdim = 2*L
  }
  X2 = cbind(D1,D2)                                # construct the matrix for the intercepts
  
  # checking that it fits
  if (FALSE) {
    eps = X1 %*% (1/model$B1) - X2 %*% c(  (model$A1 - model$r1*model$A2)/(model$B1*(1-model$r1)) ,  - (model$A4 - model$r4*model$A3)/(model$B1*(1-model$r4))   )
  }

  # construct the Z matrix of instruemts (j1,j2)
  Z = array(0,c(N,L^2))
  I = (1:N) + N*(J1-1 + L*(J2-1))
  Z[I]=1
  
  # remove empty interactions
  I = colSums(Z)!=0
  Z=Z[,I]
  Z  = Matrix(Z,sparse=T)
  
  # run LIML 
  b_liml = liml.nodep(cbind(X1,X2[,1:(2*L-1)]),Z)
  b_liml = b_liml/b_liml[norm]
  
  # extract results
  B1 = 1/b_liml[1:L]
  if (fixb==FALSE) {
    B2 = - 1/b_liml[(L+1):(2*L)]    
    A1 = - b_liml[(2*L+1):(3*L)] * B1
    A2 = c(b_liml[(3*L+1):(4*L-1)],0) * B2 
  } else {
    B2 = B1 * (1-r4)/(1-r1)
    A1 = - b_liml[(xdim):(xdim+L-1)]*B1*(1-r1)
    A2 = c(b_liml[(xdim+L):(xdim + 2*L-2)],0) * B2 *(1-r1)
  }
  
  return(list(A1=A1,A2=A2,B1=B1,B2=B2,b_liml=b_liml))
}

#' this function uses LIML to estimate a model with a given normalization
#' and we want to do it in a non-stationary way
m4.mini.linear.int <- function(Y1,Y2,J1,J2,norm=1,r1=0,r4=0,alpha=0) {
  
  L = max(J1)
  N = length(Y1)
  
  # ----- prepare the different regressors ----- #
  D1 = array(0,c(N,L))
  I = (1:N) + N*(J1-1)
  D1[I]=1
  D2 = array(0,c(N,L))
  I = (1:N) + N*(J2-1)
  D2[I]=1
  
  # because we are linear here we only need to regress on dummies
  X2 = cbind(D1[,2:L],D2) # construct the matrix for the intercepts
  fit = as.numeric(lm.fit(X2, Y1/(1-r1) - Y2/(1-r4) )$coefficients)

  A1 = c(0,fit[1:(L-1)]) * (1-r1)
  A2 = -fit[L:(2*L-1)] * (1-r4)
  B1 = rep(1-r1,L)
  B2 = rep(1-r4,L)
    
  return(list(A1=A1,A2=A2,B1=B1,B2=B2,b_liml=fit))
}





#' this function uses LIML to estimate a model with a given normalization.
#' In estimation the interaction terms B are profiled out in period 2 first, 
#' then they are used to recover the intercepts in both periods.
m4.mini.liml.int.prof <- function(Y1,Y2,J1,J2,norm=1,r1=0,r4=0) {
  
  L = max(J1)
  N = length(Y1)
  
  # ----- prepare the different regressors ----- #
  D1 = array(0,c(N,L))
  I = (1:N) + N*(J1-1)
  D1[I]=1
  D2 = array(0,c(N,L))
  I = (1:N) + N*(J2-1)
  D2[I]=1

  # construct the Z matrix of instruemts (j1,j2)
  Z = array(0,c(N,L^2))
  I = (1:N) + N*(J1-1 + L*(J2-1))
  Z[I]=1
  
  # remove empty interactions
  I = colSums(Z)!=0
  Z=Z[,I]
  Z  = Matrix(Z,sparse=T)
  
  # profiling, we project on D2Y2 and intercepts
  D2Y2 = cbind(D2*spread(Y2,2,L))
  D1Y1 = cbind(D1*spread(Y1,2,L),D1,D2[,2:L])
  X = D2Y2 - D1Y1 %*% (ginv(D1Y1) %*% D2Y2)
  X = Matrix(X,sparse=T)
  
  # Run LIML, get the Bs
  b_liml  = liml.nodep(X,Z)
  b_liml = b_liml/b_liml[norm] * (1-r1)/(1-r4) # we normalize (1-r1)*B1[norm]=1
  B2 = 1/b_liml[1:L]
  B1 = (1-r1)*B2/(1-r4)

  #b_liml  = b_liml/b_liml[norm]
  #B2 = 1/b_liml[1:L]
  #B1 = B2 # we are impossing stationarity in this estimator
   
  # finally extract the intercepts, this is a linear regression
  X = cbind(D1,-D2[,1:(L-1)])
  Y = (D1*spread(Y1,2,L)) %*% (1/B1) - (D2*spread(Y2,2,L)) %*% (1/B2)
  fit2 = lm.fit(X,Y)
  
  # testing moment condition
  if (FALSE) {  
    (D1*spread(Y1,2,L)) %*% (1/((1-r1)*model$B1)) - (D2*spread(Y2,2,L))%*% (1/((1-r4)*model$B1)) 
    (model$A1 - r1*model$A2)/(1/((1-r1)*model$B1)) - (model$A4 - r4*model$A3)/(1/((1-r4)*model$B1))
  }
  
  # --------- extract the interaction terms ---------- #
  A1 = as.numeric(fit2$coefficients[1:L]) * B1 
  A2 = c(as.numeric(fit2$coefficients[(L+1):(2*L-1)]),0) * B2
  
  if (any(is.na(A1*A2))) warnings("A1 or A2 contains NA values");
  if (any(is.na(B1*B2))) warnings("A1 or A2 contains NA values");
  
  return(list(A1=A1,A2=A2,B1=B1,B2=B2,b_liml=b_liml))
}


#' Estimates minimodel on 4 periods
#' @export
m4.mini.estimate <- function(jdata,sdata,r1,r4,model0=c(),method="ns",norm=1,alpha=0) {

  # --------- use LIML on movers to get A1,B1,A2,B2 ---------- #
  Y1=jdata$y1;Y2=jdata$y2;Y3=jdata$y3;Y4=jdata$y4;J1=jdata$j1;J2=jdata$j2;
  nf = max(J1)
  
  if (method=="ns") {
    rliml = m4.mini.liml.int(Y1-r1*Y2,Y4-r4*Y3,J1,J2,norm=norm,fixb = F,r1=r1,r4=r4,alpha=alpha)
  } else if (method=="prof") {
    rliml  = m4.mini.liml.int.prof(Y1-r1*Y2,Y4-r4*Y3,J1,J2,norm=norm,r1=r1,r4=r4)
    #rliml2 = m4.mini.liml.int(Y1-r1*Y2,Y4-r4*Y3,J1,J2,norm=norm,fixb = F,r1=r1,r4=r4,alpha=alpha)
  } else if (method=="linear") {
    rliml  = m4.mini.linear.int(Y1-r1*Y2,Y4-r4*Y3,J1,J2,norm=norm,r1=r1,r4=r4)
  } else {
    rliml = m4.mini.liml.int(Y1-r1*Y2,Y4-r4*Y3,J1,J2,norm=norm,fixb = T,r1=r1,r4=r4,alpha=alpha)
  }
  AA1 = rliml$A1 # this returns A1 - r1*A2
  AA2 = rliml$A2 # this returns A4 - r1*A3
  BB1 = rliml$B1 # this returns B1 - r1*B2
  BB2 = rliml$B2 # this returns B4 - r4*B3
  N   = length(Y1)
  Nm  = acast(jdata[,.N,list(j1,j2)],j1~j2,value.var ="N")
  Ns  = sdata[,.N,j1][order(j1)][,N]
  
  assert_notna(AA1=AA1,AA2=AA2,BB1=BB1,BB2=BB2)

  # get E[alpha|l1,l2]
  M1 = acast(jdata[,mean(y1-r1*y2),list(j1,j2)],j1~j2,fill=0,value.var ="V1")
  M2 = acast(jdata[,mean(y4-r4*y3),list(j1,j2)],j1~j2,fill=0,value.var ="V1")
  EEm = ( M1 - spread(AA1,2,nf)  )/(2*spread(BB1,2,nf)) + ( M2 - spread(AA2,1,nf)  )/(2*spread(BB2,1,nf))

  assert_notna(EEm=EEm)
  
  # get A,B,C using MEAN RESTRICTIONS
  setkey(jdata,j1,j2)
  YY1 = c(acast(jdata[,mean(y1),list(j1,j2)],j2~j1,fill=0,value.var="V1"))
  YY2 = c(acast(jdata[,mean(y2),list(j1,j2)],j2~j1,fill=0,value.var="V1"))
  YY3 = c(acast(jdata[,mean(y3),list(j1,j2)],j2~j1,fill=0,value.var="V1"))
  YY4 = c(acast(jdata[,mean(y4),list(j1,j2)],j2~j1,fill=0,value.var="V1"))
  W   = c(acast(jdata[,.N,list(j1,j2)],j2~j1,fill=0,value.var="N"))

  XX1 = array(0,c(nf^2,10*nf-2))
  XX2 = array(0,c(nf^2,10*nf-2))
  XX3 = array(0,c(nf^2,10*nf-2))
  XX4 = array(0,c(nf^2,10*nf-2))
  L = nf

  for (l1 in 1:nf) for (l2 in 1:nf) {
    ll = l2 + nf*(l1 -1)
    Ls = 0; 

    # the A
    XX1[ll,Ls+l1] = 1;Ls= Ls+L
    XX2[ll,Ls+l1] = 1;Ls= Ls+L
    XX3[ll,Ls+l2] = 1;Ls= Ls+L
    XX4[ll,Ls+l2] = 1;Ls= Ls+L

    # the B
    XX1[ll,Ls+l1] = EEm[l1,l2];Ls= Ls+L
    XX2[ll,Ls+l1] = EEm[l1,l2];Ls= Ls+L
    XX3[ll,Ls+l2] = EEm[l1,l2];Ls= Ls+L
    XX4[ll,Ls+l2] = EEm[l1,l2];Ls= Ls+L

    # the C
    if (l2>1) {
      XX1[ll,Ls+l2-1] = r1;
      XX2[ll,Ls+l2-1] = 1 ;
    }
    Ls= Ls+L-1

    if (l1>1) {
      XX3[ll,Ls+l1-1] = 1;
      XX4[ll,Ls+l1-1] = r4 ;
    }
  }

  # construct the constraints
  CM1 = cbind(diag(nf),-r1*diag(nf), array(0,c(nf,8*nf-2)))
  CM2 = cbind(array(0,c(nf,2*nf)), -r4*diag(nf),diag(nf), array(0,c(nf,6*nf-2)) )
  CM3 = cbind(array(0,c(nf,4*nf)), diag(nf),-r1*diag(nf), array(0,c(nf,4*nf-2)))
  CM4 = cbind(array(0,c(nf,6*nf)), -r4*diag(nf),diag(nf), array(0,c(nf,2*nf-2)))

  XX = rbind(XX1,XX2,XX3,XX4)
  YY = c(YY1,YY2,YY3,YY4)
  CM = rbind(CM1,CM2,CM3,CM4)
  C0 = c(AA1,AA2,BB1,BB2)
  meq = 4*nf
  
  # add 2 contraints, B1=B2 and B3=B4
  if (method %in% c("prof","fixb","linear")) {
    CM5 =  cbind(array(0,c(nf,4*nf)), diag(nf), -diag(nf), array(0,c(nf,4*nf-2)))
    CM6 =  cbind(array(0,c(nf,6*nf)), diag(nf), -diag(nf), array(0,c(nf,2*nf-2)))
    CM = rbind(CM,CM5,CM6)
    C0 = c(C0,rep(0,2*L))
    meq = meq+2*L
  }
  
  rs = lm.wfitc(XX,YY,c(W,W,W,W),CM,C0,meq)$solution
  
  Ls = 0;
  A1 = rs[(Ls+1):(Ls+L)];Ls=Ls+L
  A2 = rs[(Ls+1):(Ls+L)];Ls=Ls+L
  A3 = rs[(Ls+1):(Ls+L)];Ls=Ls+L
  A4 = rs[(Ls+1):(Ls+L)];Ls=Ls+L
  B1 = rs[(Ls+1):(Ls+L)];Ls=Ls+L
  B2 = rs[(Ls+1):(Ls+L)];Ls=Ls+L
  B3 = rs[(Ls+1):(Ls+L)];Ls=Ls+L
  B4 = rs[(Ls+1):(Ls+L)];Ls=Ls+L
  C2 = c(0,rs[(Ls+1):(Ls+L-1)]);Ls=Ls+L-1
  C3 = c(0,rs[(Ls+1):(Ls+L-1)]);Ls=Ls+L-1
  model = list(A1=A1,A2=A2,A3=A3,A4=A4,B1=B1,B2=B2,B3=B3,B4=B4,C2=C2,C3=C3,EEm=EEm,r1=r1,r4=r4,nf=nf)
  assert_notna(A1=A1,A2=A2,A3=A3,A4=A4,B1=B1,B2=B2,B3=B3,B4=B4,C2=C2,C3=C3)
  
  res_mean_stayer = m4.mini.getstayers(jdata,sdata,model)
  model$Em=res_mean_stayer$Em
  model$A2s=res_mean_stayer$A2s
  model$A3s=res_mean_stayer$A3s
  
  # get the variances
  res_var = m4.mini.extract.var(jdata,sdata,model,r1,r4)
  model$Esd    = res_var$Esd
  model$EEsd   = res_var$EEsd
  model$nu1_sd = res_var$nu1_sd
  model$nu2_sd = res_var$nu2_sd
  model$eps2m_sd = res_var$eps2m_sd
  model$eps3m_sd = res_var$eps3m_sd
  model$eps2s_sd = res_var$eps2s_sd
  model$eps3s_sd = res_var$eps3s_sd
  model$rho32m   = res_var$rho32m
  model$rho32s   = res_var$rho32s
  assert_notna(Esd=model$Esd,EEsd=model$EEsd,
               nu1_sd=model$nu1_sd,
               nu2_sd=model$nu2_sd,
               eps2s_sd=model$eps2s_sd,
               eps3s_sd=model$eps3s_sd,
               eps2m_sd=model$eps2m_sd,
               eps3m_sd=model$eps3m_sd)
  
  model$Ns=Ns
  model$Nm=Nm

  # compute vdec
  if ("k_imp" %in% names(sdata)) sdata[,k_imp:=NULL];
  sdata[,y1:=y1_bu]
  
  # simulate movers/stayers, and combine
  NNm = model$Nm; NNs = model$Ns
  stayer_share = sum(NNs)/(sum(NNs)+sum(NNm))
  model$vdec   = m4.mini.vdec(model,1e6,stayer_share,"y2")
  
  if (length(model0)>0) {
    rr = addmom(AA1,model0$A1 - model0$r1*model0$A2,"A1 - r1*A2")
    rr = addmom(AA2,model0$A4 - model0$r4*model0$A3,"A4 - r4*A3",rr)
    rr = addmom(BB1,model0$B1 - model0$r1*model0$B2,"B1 - r1*B2",rr)
    rr = addmom(BB2,model0$B4 - model0$r4*model0$B3,"B4 - r4*B3",rr)
    rr = addmom(EEm,model0$EEm,"E(alpha|l1,l2)",rr)
    rr = addmom(A1,model0$A1,"A1",rr)
    rr = addmom(A2,model0$A2,"A2",rr)
    rr = addmom(A3,model0$A3,"A3",rr)
    rr = addmom(A4,model0$A4,"A4",rr)
    rr = addmom(B1,model0$B1,"B1",rr)
    rr = addmom(B2,model0$B2,"B2",rr)
    rr = addmom(B3,model0$B3,"B3",rr)
    rr = addmom(B4,model0$B4,"B4",rr)
    rr = addmom(C2,model0$C2,"C2",rr)
    rr = addmom(C3,model0$C3,"C3",rr)
    rr = addmom(model$Em,model0$Em,"Em",rr)
    rr = addmom(model$Esd,model0$Esd,"Esd",rr)
    rr = addmom(model$EEsd,model0$EEsd,"EEsd",rr)
    rr = addmom(model$nu1_sd,model0$nu1_sd,"nu1_sd",rr)
    
    #rr = addmom(model$nu2_sd,model0$nu2_sd,"nu2_sd",rr)
    print(ggplot(rr,aes(x=val2,y=val1,color=type)) + geom_point() + facet_wrap(~name,scale="free") + theme_bw() + geom_abline(linetype=2))
    catf("cor with true model:%f",cor(rr$val1,rr$val2))
  }

  return(model)
}

#' Dynamic Linear model
#'
#' @param jdata 
#' @param sdata 
#' @param r1 
#' @param r4 
#' @param model0 
#' @param norm which firm to normalize to 1
#'
#' @export
m4.minilin.estimate <- function(jdata,sdata,r1,r4,model0=c(),norm=1) {
  
  # --------- use LIML on movers to get A1,B1,A2,B2 ---------- #
  Y1=jdata$y1;Y2=jdata$y2;Y3=jdata$y3;Y4=jdata$y4;J1=jdata$j1;J2=jdata$j2;
  nf = max(J1)
  N   = length(Y1)

  # get A,B,C using MEAN RESTRICTIONS
  setkey(jdata,j1,j2)
  YY1 = as.numeric(acast(jdata[,mean(y1),list(j1,j2)],j2~j1,value.var = "V1",fill=0))
  YY2 = as.numeric(acast(jdata[,mean(y2),list(j1,j2)],j2~j1,value.var = "V1",fill=0))
  YY3 = as.numeric(acast(jdata[,mean(y3),list(j1,j2)],j2~j1,value.var = "V1",fill=0))
  YY4 = as.numeric(acast(jdata[,mean(y3),list(j1,j2)],j2~j1,value.var = "V1",fill=0))
  W   = as.numeric(acast(jdata[,list(V1=.N),list(j1,j2)],j2~j1,value.var = "V1",fill=0))
  Nm  = acast(jdata[,.N,list(j1,j2)],j1~j2,value.var ="N",fill=0)
  Ns  = sdata[,.N,j1][order(j1)][,N]

  XX1 = array(0,c(nf^2,6*nf-2 + nf^2))
  XX2 = array(0,c(nf^2,6*nf-2 + nf^2))
  XX3 = array(0,c(nf^2,6*nf-2 + nf^2))
  XX4 = array(0,c(nf^2,6*nf-2 + nf^2))
  L = nf

  for (l1 in 1:nf) for (l2 in 1:nf) {
    ll = l2 + nf*(l1 -1)
    Ls = 0; 

    # the A
    XX1[ll,Ls+l1] = 1;Ls= Ls+L
    XX2[ll,Ls+l1] = 1;Ls= Ls+L
    XX3[ll,Ls+l2] = 1;Ls= Ls+L
    XX4[ll,Ls+l2] = 1;Ls= Ls+L

    # the EEm
    XX1[ll,Ls+ll] = 1;
    XX2[ll,Ls+ll] = 1;
    XX3[ll,Ls+ll] = 1;
    XX4[ll,Ls+ll] = 1;Ls= Ls+L^2

    # the C
    if (l2>1) {
      XX1[ll,Ls+l2-1] = r1;
      XX2[ll,Ls+l2-1] = 1 ;
    }
    Ls= Ls+L-1

    if (l1>1) {
      XX3[ll,Ls+l1-1] = 1;
      XX4[ll,Ls+l1-1] = r4 ;
    }
  }
    
  XX = rbind(XX1,XX2,XX3,XX4)
  XX = XX[,2:ncol(XX)] # removing the first intercept
  YY = c(YY1,YY2,YY3,YY4)
  
  rs = coef(lm.wfit(XX,YY,c(W,W,W,W)))
  
  Ls = 0;
  A1 = c(0,as.numeric(rs[(Ls+1):(Ls+L-1)]));Ls=Ls+L-1
  A2 = as.numeric(rs[(Ls+1):(Ls+L)]);Ls=Ls+L
  A3 = as.numeric(rs[(Ls+1):(Ls+L)]);Ls=Ls+L
  A4 = as.numeric(rs[(Ls+1):(Ls+L)]);Ls=Ls+L
  EEm= t(array( as.numeric(rs[(Ls+1):(Ls+L^2)]),c(L,L)));Ls=Ls+L^2
  C2 = c(0,as.numeric(rs[(Ls+1):(Ls+L-1)]));Ls=Ls+L-1
  C3 = c(0,as.numeric(rs[(Ls+1):(Ls+L-1)]));Ls=Ls+L-1
  
  model = list(A1=A1,A2=A2,A3=A3,A4=A4,B1=rep(1,nf),B2=rep(1,nf),B3=rep(1,nf),B4=rep(1,nf),C2=C2,C3=C3,EEm=EEm,r1=r1,r4=r4,nf=nf)
  res_mean_stayer = m4.mini.getstayers(jdata,sdata,model)
  model$Em=res_mean_stayer$Em
  model$A2s=res_mean_stayer$A2s
  model$A3s=res_mean_stayer$A3s
  
  # get the variances
  res_var = m4.mini.getvar.movers.unc(jdata,r1,r4)
  model$Esd=res_var$Esd
  model$EEsd=res_var$EEsd
  model$nu1_sd=res_var$nu1_sd
  model$nu2_sd=res_var$nu2_sd
  model$fit1 = res_var$fit1
  model$fit2 = res_var$fit2
  model$Evar=res_var$Evar
  
  # get the variances from the full stayers variance structure
  model$fit3 = m4.mini.getvar.stayers(jdata,sdata,model,r1,r4)
  
  model$Ns=Ns
  model$Nm=Nm
  
  model = model.mini.rho32.stayers(jdata,sdata,model)
  model = model.mini.rho32.movers(jdata,sdata,model)
  
  if (length(model0)>0) {
    rr = addmom(A1,model0$A1,"A1")
    rr = addmom(A2,model0$A2,"A2",rr)
    rr = addmom(A3,model0$A3,"A3",rr)
    rr = addmom(A4,model0$A4,"A4",rr)
    rr = addmom(C2,model0$C2,"C2",rr)
    rr = addmom(C3,model0$C3,"C3",rr)
    rr = addmom(EEm,model0$EEm,"EEm",rr)
    rr = addmom(model$Em,model0$Em,"Em",rr)
    rr = addmom(model$A2s,model0$A2s,"A2s",rr)
    rr = addmom(model$A3s,model0$A3s,"A3s",rr)
    rr = addmom(model$Esd,model0$Esd,"Esd",rr)
    rr = addmom(model$EEsd,model0$EEsd,"EEsd",rr)
    rr = addmom(model$nu1_sd,model0$nu1_sd,"nu1_sd",rr)
    rr = addmom(model$nu2_sd,model0$nu2_sd,"nu2_sd",rr)
    print(ggplot(rr,aes(x=val2,y=val1,color=type)) + geom_point() + facet_wrap(~name,scale="free") + theme_bw() + geom_abline(linetype=2))
    catf("cor with true model:%f",cor(rr$val1,rr$val2))
  }
  
  return(model)
}

#' This function implements quadratic programing described in the appendix
#' to recover nu_1, nu_2, eps_sd_1, eps_sd_2 and the rho32m and rho32m
#' this uses the B parameters from the minimodel.
m4.mini.extract.var <- function(jdata,sdata,model,r1,r4) {

  # -------- USE MOVERS ---------- #
  setkey(jdata,j1,j2)
  YY0_11 = c(acast(jdata[, mvar(y1-r1*y2)           ,list(j1,j2)],j2~j1,fill=0,value.var="V1"))
  YY1_22 = c(acast(jdata[, mvar(y2)                 ,list(j1,j2)],j2~j1,fill=0,value.var="V1"))
  YY2_33 = c(acast(jdata[, mvar(y3)                 ,list(j1,j2)],j2~j1,fill=0,value.var="V1"))
  YY3_44 = c(acast(jdata[, mvar(y4-r4*y3)           ,list(j1,j2)],j2~j1,fill=0,value.var="V1"))
  YY4_12 = c(acast(jdata[, mcov(y1-r1*y2, y2)       ,list(j1,j2)],j2~j1,fill=0,value.var="V1"))
  YY5_13 = c(acast(jdata[, mcov(y1-r1*y2, y3)       ,list(j1,j2)],j2~j1,fill=0,value.var="V1"))
  YY6_14 = c(acast(jdata[, mcov(y1-r1*y2, y4-r4*y3) ,list(j1,j2)],j2~j1,fill=0,value.var="V1"))
  YY8_24 = c(acast(jdata[, mcov(y2      , y4-r4*y3) ,list(j1,j2)],j2~j1,fill=0,value.var="V1"))
  YY9_34 = c(acast(jdata[, mcov(y3      , y4-r4*y3) ,list(j1,j2)],j2~j1,fill=0,value.var="V1"))
  W      = c(acast(jdata[,.N,list(j1,j2)],j2~j1,fill=0,value.var="N"))

  nf = model$nf
  XX0_11 = array(0,c(nf^2,3*nf^2 + 2*nf))
  XX1_22 = array(0,c(nf^2,3*nf^2 + 2*nf))
  XX2_33 = array(0,c(nf^2,3*nf^2 + 2*nf))
  XX3_44 = array(0,c(nf^2,3*nf^2 + 2*nf))
  XX4_12 = array(0,c(nf^2,3*nf^2 + 2*nf))
  XX5_13 = array(0,c(nf^2,3*nf^2 + 2*nf))
  XX6_14 = array(0,c(nf^2,3*nf^2 + 2*nf))
  #XX7_23 = array(0,c(nf^2,4*nf^2 + 2*nf))
  XX8_24 = array(0,c(nf^2,3*nf^2 + 2*nf))
  XX9_34 = array(0,c(nf^2,3*nf^2 + 2*nf))
  L = nf

  B1 = model$B1
  B2 = model$B2
  B3 = model$B3
  B4 = model$B4

  BB1 = B1 - model$r1*B2
  BB2 = B4 - model$r4*B3

  for (l1 in 1:nf) for (l2 in 1:nf) {
    ll = l2 + nf*(l1 -1)
    Ls = 0; 

    if (W[ll]>0) {
      # the Var(alpha|l1,l2)
      # on diagonal
      XX0_11[ll,ll] = BB1[l1]^2
      XX1_22[ll,ll] = B2[l1]^2
      XX2_33[ll,ll] = B3[l2]^2
      XX3_44[ll,ll] = BB2[l2]^2
      # off diagonal
      XX4_12[ll,ll] = BB1[l1]*B2[l1]
      XX5_13[ll,ll] = BB1[l1]*B3[l2]
      XX6_14[ll,ll] = BB1[l1]*BB2[l2]
      XX8_24[ll,ll] = B2[l1] *BB2[l2]
      XX9_34[ll,ll] = B3[l2] *BB2[l2]
      Ls= Ls+L^2
  
      #  var(epsilon) 
      XX1_22[ll,Ls+ll] = 1;Ls= Ls+L^2
      XX2_33[ll,Ls+ll] = 1;Ls= Ls+L^2
  
      # the var(nu | l) 
      XX0_11[ll,Ls+l1] = 1;Ls= Ls+L
      XX3_44[ll,Ls+l2] = 1;Ls= Ls+L
    } else {
      
      # constraint the alpha to be zero 
      XX0_11[ll,ll] = 1
      XX1_22[ll,ll] = 1
      XX2_33[ll,ll] = 1
      XX3_44[ll,ll] = 1
      Ls= Ls+L^2
      XX1_22[ll,Ls+ll] = 1;Ls= Ls+L^2
      XX2_33[ll,Ls+ll] = 1;Ls= Ls+L^2
      W[ll] = 0.1
    }
  }

  XX = rbind(XX0_11,XX1_22,XX2_33,XX3_44,XX4_12,XX5_13,XX6_14,XX8_24,XX9_34)
  YY =     c(YY0_11,YY1_22,YY2_33,YY3_44,YY4_12,YY5_13,YY6_14,YY8_24,YY9_34)

  fit1 = lm.wfitnn(XX,YY,c(W,W,W,W,W,W,W,W,W))  
  res  = fit1$solution
  obj1 = t( YY - XX %*% res ) %*% diag(c(W,W,W,W,W,W,W,W,W)) %*% ( YY - XX %*% res )
 
  Ls = 0
  EEsd     = t(sqrt(pmax(array((res)[(Ls+1):(Ls+nf^2)],c(nf,nf)),0))); Ls = Ls + nf^2
  eps2m_sd = t(sqrt(pmax(array((res)[(Ls+1):(Ls+nf^2)],c(nf,nf)),0))); Ls = Ls + nf^2
  eps3m_sd = t(sqrt(pmax(array((res)[(Ls+1):(Ls+nf^2)],c(nf,nf)),0))); Ls = Ls + nf^2
  #cov23_sd = t(sqrt(pmax(array((res)[(Ls+1):(Ls+nf^2)],c(nf,nf)),0))); Ls = Ls + nf^2
  nu1_sd   = sqrt(pmax((res)[(Ls+1):(Ls+nf)],0)); Ls = Ls +nf
  nu2_sd   = sqrt(pmax((res)[(Ls+1):(Ls+nf)],0)); Ls = Ls +nf

  # get rho32m
  fit_rho = jdata[, list(y=cov(y2 , y3) - B2[j1] *B3[j2] * EEsd[j1,j2]^2 ,x= eps2m_sd[j1,j2] * eps3m_sd[j1,j2]),list(j1,j2)][, lm(y~0+x)]
  rho32m = coef(fit_rho)[[1]]

  # -------------  USE STAYERS ---------------- #
  setkey(sdata,j1)
  YY0_11 = sdata[, var(y1-r1*y2)           ,list(j1)][,V1] - nu1_sd^2
  YY1_22 = sdata[, var(y2)                 ,list(j1)][,V1]
  YY2_33 = sdata[, var(y3)                 ,list(j1)][,V1]
  YY3_44 = sdata[, var(y4-r4*y3)           ,list(j1)][,V1] - nu2_sd^2
  YY4_12 = sdata[, cov(y1-r1*y2, y2)       ,list(j1)][,V1]
  YY5_13 = sdata[, cov(y1-r1*y2, y3)       ,list(j1)][,V1]
  YY6_14 = sdata[, cov(y1-r1*y2, y4-r4*y3) ,list(j1)][,V1]
  #YY7_23 = sdata[, cov(y2      , y3)       ,list(j1)][,V1]
  YY8_24 = sdata[, cov(y2      , y4-r4*y3) ,list(j1)][,V1]
  YY9_34 = sdata[, cov(y3      , y4-r4*y3) ,list(j1)][,V1]
  W   = sdata[,.N,list(j1)][,N]

  nf = model$nf
  XX0_11 = array(0,c(nf,3*nf))
  XX1_22 = array(0,c(nf,3*nf))
  XX2_33 = array(0,c(nf,3*nf))
  XX3_44 = array(0,c(nf,3*nf))
  XX4_12 = array(0,c(nf,3*nf))
  XX5_13 = array(0,c(nf,3*nf))
  XX6_14 = array(0,c(nf,3*nf))
  #XX7_23 = array(0,c(nf,4*nf))
  XX8_24 = array(0,c(nf,3*nf))
  XX9_34 = array(0,c(nf,3*nf))
  L = nf

  B1 = model$B1
  B2 = model$B2
  B3 = model$B3
  B4 = model$B4

  BB1 = B1 - model$r1*B2
  BB2 = B4 - model$r4*B3

  for (l1 in 1:nf) {
    Ls = 0; 

    # the Var(alpha|l1)
    # on diagonal
    XX0_11[l1,l1] = BB1[l1]^2
    XX1_22[l1,l1] = B2[l1]^2
    XX2_33[l1,l1] = B3[l1]^2
    XX3_44[l1,l1] = BB2[l1]^2
    # off diagonal
    XX4_12[l1,l1] = BB1[l1]*B2[l1]
    XX5_13[l1,l1] = BB1[l1]*B3[l1]
    XX6_14[l1,l1] = BB1[l1]*BB2[l1]
    #XX7_23[l1,l1] = B2[l1] *B3[l1]
    XX8_24[l1,l1] = B2[l1] *BB2[l1]
    XX9_34[l1,l1] = B3[l1] *BB2[l1]
    Ls= Ls+L

    #  var(epsilon), cov by l1
    XX1_22[l1,Ls+l1] = 1;Ls= Ls+L
    XX2_33[l1,Ls+l1] = 1;Ls= Ls+L
    #XX7_23[l1,Ls+l1] = 1;Ls= Ls+L
  }

  #XX = rbind(XX0_11,XX1_22,XX2_33,XX3_44,XX4_12,XX5_13,XX6_14,XX7_23,XX8_24,XX9_34)
  #YY =     c(YY0_11,YY1_22,YY2_33,YY3_44,YY4_12,YY5_13,YY6_14,YY7_23,YY8_24,YY9_34)
  XX = rbind(XX0_11,XX1_22,XX2_33,XX3_44,XX4_12,XX5_13,XX6_14,XX8_24,XX9_34)
  YY =     c(YY0_11,YY1_22,YY2_33,YY3_44,YY4_12,YY5_13,YY6_14,YY8_24,YY9_34)

  fit2 = lm.wfitnn(XX,YY,c(W,W,W,W,W,W,W,W,W))  
  res  = fit2$solution
  obj2 = t( YY - XX %*% res ) %*% diag(c(W,W,W,W,W,W,W,W,W)) %*% ( YY - XX %*% res )
 
  Ls = 0
  Esd       = t(sqrt(pmax((res)[(Ls+1):(Ls+nf)],0))); Ls = Ls + L
  eps2s_sd  = t(sqrt(pmax((res)[(Ls+1):(Ls+nf)],0))); Ls = Ls + L
  eps3s_sd  = t(sqrt(pmax((res)[(Ls+1):(Ls+nf)],0))); Ls = Ls + L

  # get rho32s
  fit_rho = sdata[, list(y=cov(y2 , y3) - B2[j1] *B3[j1] * Esd[j1]^2 ,x= eps2s_sd[j1] * eps3s_sd[j1]),list(j1)][, lm(y~0+x)]
  rho32s = coef(fit_rho)[[1]]
  
  return(list(Esd=Esd,EEsd=EEsd,nu1_sd=nu1_sd,nu2_sd=nu2_sd,eps2m_sd=eps2m_sd,eps3m_sd=eps3m_sd,fit1=obj1,fit2=obj2,
              eps2s_sd=eps2s_sd,eps3s_sd=eps3s_sd,rho32s=rho32s,rho32m=rho32m,Evar=res))
}

#' Computes the variance decomposition by simulation
#' @export
m4.mini.vdec <- function(model,nsim,stayer_share=1,ydep="y2") {
  
  # simulate movers/stayers, and combine
  NNm = model$Nm
  NNs = model$Ns
  NNm[!is.finite(NNm)]=0
  NNs[!is.finite(NNs)]=0
  
  NNs = round(NNs*nsim*stayer_share/sum(NNs))
  NNm = round(NNm*nsim*(1-stayer_share)/sum(NNm))
  flog.info("computing var decomposition with ns=%i nm=%i",sum(NNs),sum(NNm))
  
  # we simulate from the model both movers and stayers
  sdata.sim = m4.mini.simulate.stayers(model,NNs)
  jdata.sim = m4.mini.simulate.movers(model,NNm)
  sdata.sim = rbind(sdata.sim[,list(j1,k=alpha,y1,y2,y3,y4)],jdata.sim[,list(j1,k=alpha,y1,y2,y3,y4)])
  proj_unc  = lin.proja(sdata.sim,ydep,"k","j1");  
  
  return(proj_unc)
}



test.m4.mini.getvar <- function() {
  model = m4.mini.new(10,r1=0.4,r4=0.7,eps_sd_factor = 1)
  NNm   = array(20000/(model$nf^2),c(model$nf,model$nf))
  NNs  = array(100000/model$nf,model$nf)
  jdata = m4.mini.simulate.movers(model,NNm)
  sdata = m4.mini.simulate.stayers(model,NNs)
  
  model1 = m4.mini.getvar(jdata,sdata,model,r1=model$r1,r4=model$r4)
  plot(model$nu1_sd,model1$nu1_sd)
  plot(model$nu2_sd,model1$nu2_sd)
}




model.mini.rho32.movers <- function(jdatae,sdata,model) {
  # get the variances of the epsilons
  eps2_sd = sqrt(pmax(acast(jdatae[, var(y2), list(j1,j2)],j1~j2,value.var = "V1") - model$EEsd^2,0))
  eps3_sd = sqrt(pmax(acast(jdatae[, var(y3), list(j1,j2)],j1~j2,value.var = "V1") - model$EEsd^2,0))
  
  # get the correlation  
  rtmp = jdatae[,list(y = cov(y2,y3) - model$EEsd[j1,j2]^2, x = eps2_sd[j1,j2]*eps3_sd[j1,j2],.N),list(j1,j2)]
  fit = lm(y~0+x,rtmp,weights = rtmp$N)

  model$rho32m   = coef(fit)[[1]]
  model$eps2m_sd = eps2_sd
  model$eps3m_sd = eps3_sd
  
  return(model)
}

model.mini.rho32.stayers <- function(jdatae,sdata,model) {
  # get the variances of the epsilons
  eps2s_sd = sqrt(pmax(sdata[, var(y2), j1][,V1] - model$Esd^2,0))
  eps3s_sd = sqrt(pmax(sdata[, var(y3), j1][,V1] - model$Esd^2,0))
  
  # get the correlation
  rtmp = sdata[,list(y = cov(y2,y3) - model$Esd[j1]^2, x = eps2s_sd[j1]*eps3s_sd[j1],.N),j1]
  fit = lm(y~0+x,rtmp,weights = rtmp$N)
  
  model$rho32s   = coef(fit)[[1]]
  model$eps3s_sd = eps3s_sd
  model$eps2s_sd = eps2s_sd
  
  return(model)
}


# extract Var(alpha|l1,l2) and Var(epsion|l)
m4.mini.getvar.stayers <- function(jdata,sdata,model,r1,r4) {

  B1 = model$B1
  B2 = model$B2
  B3 = model$B3
  B4 = model$B4

  BB1 = B1 - r1*B2
  BB2 = B4 - r4*B3
  
  Esd = model$Esd
  nu1_sd = model$nu1_sd
  nu2_sd = model$nu2_sd
  
  # USE MOVERS TO GET Var(epsilon|l1)
  setkey(sdata,j1)
  R1  = sdata[, var(y1-r1*y2)            - BB1[j1]^2       * Esd[j1]^2 - nu1_sd[j1]^2 ,list(j1)][,V1]
  R2  = sdata[, cov(y1-r1*y2, y2)        - BB1[j1]*B2[j1]  * Esd[j1]^2                ,list(j1)][,V1]
  R3  = sdata[, cov(y1-r1*y2, y3)        - BB1[j1]*B3[j1]  * Esd[j1]^2                ,list(j1)][,V1]
  R4  = sdata[, cov(y1-r1*y2, y4-r4*y3)  - BB1[j1]*BB2[j1] * Esd[j1]^2                ,list(j1)][,V1]
  R5  = sdata[, cov(y2      , y4-r4*y3)  - BB2[j1]*B2[j1]  * Esd[j1]^2                ,list(j1)][,V1]
  R6  = sdata[, cov(y3      , y4-r4*y3)  - BB2[j1]*B3[j1]  * Esd[j1]^2                ,list(j1)][,V1]
  R7  = sdata[, var(y4-r4*y3)            - BB2[j1]^2       * Esd[j1]^2 - nu2_sd[j1]^2 ,list(j1)][,V1]
  W   = sdata[,.N,list(j1)][,N]

  obj1 = sum( c(R1,R2,R3,R4,R5,R6,R7)^2 * c(W,W,W,W,W,W,W) )
  
  return(obj1)
}



#' we extract E(alpha|l) and A2s and A3s jointly
m4.mini.getstayers <- function(jdata,sdata,model) {

  B1 = model$B1
  B2 = model$B2
  B3 = model$B3
  B4 = model$B4
  A1 = model$A1
  A2 = model$A2
  A3 = model$A3
  A4 = model$A4
  
  BB1 = B1 - model$r1*B2
  BB2 = B4 - model$r4*B3
  AA1 = A1 - model$r1*A2
  AA2 = A4 - model$r4*A3
  r1  = model$r1
  r4  = model$r4
  nf=model$nf
  
  setkey(sdata,j1)
  YY1 = sdata[, mean(y1 -r1*y2 -AA1[j1])     ,list(j1)][,V1]
  YY2 = sdata[, mean(y2)                     ,list(j1)][,V1]
  YY3 = sdata[, mean(y3)                     ,list(j1)][,V1]
  YY4 = sdata[, mean(y4-r4*y3 - AA2[j1])     ,list(j1)][,V1]
  W   = sdata[,.N,list(j1)][,N]
  
  XX1 = array(0,c(nf,3*nf))
  XX2 = array(0,c(nf,3*nf))
  XX3 = array(0,c(nf,3*nf))
  XX4 = array(0,c(nf,3*nf))
  L = nf
  
  for (l1 in 1:nf)  {
    XX1[l1,l1]   = BB1[l1]
    XX2[l1,l1]   = B2[l1]
    XX3[l1,l1]   = B3[l1]
    XX4[l1,l1]   = BB2[l1]
    
    # A2s and A3s
    XX2[l1,L+l1]   = 1
    XX3[l1,L+L+l1] = 1
  }
  
  XX = rbind(XX1,XX2,XX3,XX4)
  YY =     c(YY1,YY2,YY3,YY4)

  fit2 = lm.wfit(XX,YY,c(W,W,W,W))
  rs = coef(fit2)
  
  Ls = 0;
  Em  = as.numeric(rs[(Ls+1):(Ls+L)]);Ls=Ls+L
  A2s = as.numeric(rs[(Ls+1):(Ls+L)]);Ls=Ls+L
  A3s = as.numeric(rs[(Ls+1):(Ls+L)]);Ls=Ls+L
  
  return(list(Em=Em,A2s=A2s,A3s=A3s))
}


m4.mini.vardec <- function(model,verbose=F) {
  # compute the common b
  B0 = wtd.mean(model$B2*model$Esd^2,model$Ns)/wtd.mean(model$Esd^2,model$Ns)
  At = model$A2s + (model$B2 - B0)*model$Em
  
  # compute the variances
  v_psi    = wtd.var(At,model$Ns) 
  v_alpha  = B0^2 * (  wtd.var(model$Em,model$Ns)  +  wtd.mean(model$Esd^2,model$Ns)   )
  v_2cov   = 2*cov.wt(cbind( At , B0 *model$Em),as.numeric(model$Ns))$cov[1,2]  
  v_eps    = wtd.mean(model$eps2s_sd^2,model$Ns)
  v_cor    = v_2cov/(2 * sqrt(v_alpha*v_psi)) 
  
  if (verbose) {
    catf("v_psi= %4.4f   v_alpha= %4.4f   2cov= %4.4f   v_eps= %4.4f \n",v_psi,v_alpha,v_2cov,v_eps)
    catf("  psi= %4.4f   v_alpha= %4.4f   2cov= %4.4f   cor  = %4.4f \n",100*v_psi/(v_psi+v_alpha+v_2cov),100*v_alpha/(v_psi+v_alpha+v_2cov),100*v_2cov/(v_psi+v_alpha+v_2cov),v_cor)
    catf("V(E(alpha|l))= %f \n", wtd.var(model$Em,model$Ns))
  }
  
  return(list(v_psi=v_psi,v_alpha=v_alpha,v_2cov=v_2cov,v_eps=v_eps))
}

m4.mini.plot <- function(m) {

  dd = data.frame()
  L = m$nf
  dd = rbind(dd,data.frame(l=1:L,y=m$A1,name="a",t=1,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$A2,name="a",t=2,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$A3,name="a",t=3,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$A4,name="a",t=4,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$A2s,name="as",t=2,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$A3s,name="as",t=3,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$C2,name="xi",t=2,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$C3,name="xi",t=3,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$B1,name="B",t=1,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$B2,name="B",t=2,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$B3,name="B",t=3,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$B4,name="B",t=4,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$Em - m$Em[1],name="alpha_l",t=2,rho1=m$r1,rho2=m$r4))

  ggplot(dd,aes(x=factor(l),y=y,group=t,color=factor(t))) + geom_line() + geom_point() +
    theme_bw() + facet_wrap(~name,scale="free")
}

#' plotting mini4 model 
#' @export 
m4.mini.plotw <- function(model,qt=6,getvals=F) {
  
  rr = data.table(l=1:length(model$A1),Em = model$Em,Esd = as.numeric(model$Esd),
                  N = model$Ns,A1=model$A1,B1=model$B1,A2s=model$A2s,B2=model$B2,B3=model$B3,B4=model$B4,A3s=model$A3s,A4=model$A4)
  
  alpha_m  = rr[, wtd.mean(Em,N)]
  alpha_sd = sqrt(rr[, wtd.mean(Esd^2,N) + wtd.var(Em,N) ])
  
  qts = qnorm( (1:qt)/(qt+1))

  rr2 = rr[, list( y1= (qts*alpha_sd + alpha_m)* B1+ A1 ,
                   y2= (qts*alpha_sd + alpha_m)* B2+ A2s,
                   y3= (qts*alpha_sd + alpha_m)* B3+ A3s,
                   y4= (qts*alpha_sd + alpha_m)* B4+ A4,
                   k=1:qt),l ]
  rr2 = melt(rr2,id.vars = c('l','k'))
  if (getvals==TRUE) {
    return(data.table(rr2))
  }
  
  ggplot(rr2,aes(x=l,y=value,color=factor(k))) + geom_line() + geom_point() + theme_bw() + facet_wrap(~variable,nrow = 2,scales="free")
}


# extract Var(alpha|l1,l2) and Var(epsion|l)
m4.mini.getvar.stayers.unc <- function(sdata,r1,r4,weights=rep(1,7)) {

  # USE MOVERS TO GET Var(epsilon|l1)
  setkey(sdata,j1)
  YY1 = sdata[, var(y1-r1*y2)           ,j1][,V1]
  YY2 = sdata[, cov(y1-r1*y2, y2)       ,j1][,V1]
  YY3 = sdata[, cov(y1-r1*y2, y3)       ,j1][,V1]
  YY4 = sdata[, cov(y1-r1*y2, y4-r4*y3) ,j1][,V1]
  YY5 = sdata[, cov(y2      , y4-r4*y3) ,j1][,V1]
  YY6 = sdata[, cov(y3      , y4-r4*y3) ,j1][,V1]
  YY7 = sdata[, var(y4-r4*y3)           ,j1][,V1]
  W   = sdata[,.N,j1][,N]

  nf  = max(sdata$j1)
  XX1 = array(0,c(nf,3*nf))
  XX2 = array(0,c(nf,3*nf))
  XX3 = array(0,c(nf,3*nf))
  XX4 = array(0,c(nf,3*nf))
  XX5 = array(0,c(nf,3*nf))
  XX6 = array(0,c(nf,3*nf))
  XX7 = array(0,c(nf,3*nf))
  L = nf

  for (l1 in 1:nf) {
    ll = l1
    Ls = 0; 

    # the Var(alpha|l1,l2)
    XX1[ll,ll] = (1-r1)^2
    XX2[ll,ll] = (1-r1)
    XX3[ll,ll] = (1-r1)
    XX4[ll,ll] = (1-r1)*(1-r4)
    XX5[ll,ll] = (1-r4)
    XX6[ll,ll] = (1-r4)
    XX7[ll,ll] = (1-r4)^2
    Ls= Ls+L

    # the var(nu|l)
    XX1[ll,Ls+l1] = 1;Ls= Ls+L
    XX7[ll,Ls+l1] = 1
  }

  XX = rbind(XX1,XX2,XX3,XX4,XX5,XX6,XX7)
  YY =     c(YY1,YY2,YY3,YY4,YY5,YY6,YY7)

  fit1 = lm.wfitnn(XX,YY,c(W*weights[1],W*weights[2],W*weights[3],W*weights[4],W*weights[5],W*weights[6],W*weights[7]))  
  res  = fit1$solution
  obj1 = t( YY - XX %*% res ) %*% diag(c(W,W,W,W,W,W,W)) %*% ( YY - XX %*% res )
 
  EEsd    = t(sqrt(pmax(array((res)[1:nf^2],c(nf,nf)),0)))
  nu1_sd = sqrt(pmax((res)[(nf^2+1):(nf^2+nf)],0))
  nu2_sd = sqrt(pmax((res)[(nf^2+nf+1):(nf^2+2*nf)],0))
  
  return(list(EEsd=EEsd,nu1_sd=nu1_sd,nu2_sd=nu2_sd,fit1=obj1,res=colMeans(rdim(YY - XX %*% res,10,7)^2)))
}

#' extract the rho parameters from stayers using variance/co-variance 
#' restrictions
#' @export
m4.mini.getvar.stayers.unc2.opt <- function(sdata,weights=c(),start=c(0.5,0.5,0.5),diff=FALSE) {
  # simple starting value strategy
  ff <- function(rhos) m4.mini.getvar.stayers.unc2(sdata,rhos[1],rhos[2],rhos[3],as.numeric(weights),diff=diff)$fit1
  rr = optim(start,ff,control=list(trace=1,REPORT=1),method="BFGS")
  res = m4.mini.getvar.stayers.unc2(sdata,rr$par[1],rr$par[2],rr$par[3],as.numeric(weights),diff=diff)
  flog.info("r1=%f r4=%f",res$r1,res$r4)
  return(res)
}





# extract Var(alpha|l1,l2) and Var(epsion|l)
m4.mini.getvar.stayers.unc2.bis <- function(sdata,r1,r4,rt,weights=c(),model0=c()) {

  # USE MOVERS TO GET Var(epsilon|l1)
  setkey(sdata,j1)
  YY1 = sdata[, var(y1)                 ,j1][,V1]
  YY2 = sdata[, var(y2)                 ,j1][,V1]
  YY3 = sdata[, var(y3)                 ,j1][,V1]
  YY4 = sdata[, var(y4)                 ,j1][,V1]
  YY5 = sdata[, cov(y1,y2)              ,j1][,V1]
  YY6 = sdata[, cov(y1,y3)              ,j1][,V1]
  YY7 = sdata[, cov(y1,y4)              ,j1][,V1]
  YY8 = sdata[, cov(y2,y3)              ,j1][,V1]
  YY9 = sdata[, cov(y2,y4)              ,j1][,V1]
  YY0 = sdata[, cov(y3,y4)              ,j1][,V1]
  W   = sdata[,.N,j1][,N]

  nf  = max(sdata$j1)
  XX1 = array(0,c(nf,5*nf))
  XX2 = array(0,c(nf,5*nf))
  XX3 = array(0,c(nf,5*nf))
  XX4 = array(0,c(nf,5*nf))
  XX5 = array(0,c(nf,5*nf))
  XX6 = array(0,c(nf,5*nf))
  XX7 = array(0,c(nf,5*nf))
  XX8 = array(0,c(nf,5*nf))
  XX9 = array(0,c(nf,5*nf))
  XX0 = array(0,c(nf,5*nf))
  L = nf

  for (l1 in 1:nf) {
    ll = l1
    Ls = 0; 

    # the Var(alpha|l1,l2)
    XX1[ll,ll] = 1
    XX2[ll,ll] = 1
    XX3[ll,ll] = 1
    XX4[ll,ll] = 1
    XX5[ll,ll] = 1
    XX6[ll,ll] = 1
    XX7[ll,ll] = 1
    XX8[ll,ll] = 1
    XX9[ll,ll] = 1
    XX0[ll,ll] = 1
    Ls= Ls+L

    # Var(nu_1|l)
    XX1[ll,Ls+l1] = 1;
    Ls= Ls+L

    # Var(nu_2|l)
    XX1[ll,Ls+l1] = r1^2;
    XX2[ll,Ls+l1] = 1;
    XX3[ll,Ls+l1] = rt^2
    XX4[ll,Ls+l1] = r1^2 * rt^2
    XX5[ll,Ls+l1] = r1
    XX6[ll,Ls+l1] = r1*rt
    XX7[ll,Ls+l1] = r1*rt*r4
    XX8[ll,Ls+l1] = rt
    XX9[ll,Ls+l1] = rt*r4
    XX0[ll,Ls+l1] = rt^2*r4 
    Ls= Ls+L

    # Var(nu_3|l)
    XX3[ll,Ls+l1] = 1
    XX4[ll,Ls+l1] = r4^2
    XX0[ll,Ls+l1] = r4
    Ls= Ls+L

    # Var(nu_4|l)
    XX4[ll,Ls+l1] = 1
  }

  XX = rbind(XX1,XX2,XX3,XX4,XX5,XX6,XX7,XX8,XX9,XX0)
  YY =     c(YY1,YY2,YY3,YY4,YY5,YY6,YY7,YY8,YY9,YY0)

  if (length(weights)==0) {
    weights =  c(W,W,W,W,W,W,W,W,W,W)
  } 
  fit1 = lm.wfitnn(XX,YY,weights)  
  res  = fit1$solution
  obj1 = t( YY - XX %*% res ) %*% diag(weights) %*% ( YY - XX %*% res )
 
  BEsd   = sqrt(pmax(res[(0*nf+1):(1*nf)],0))
  nu1_sd = sqrt(pmax(res[(1*nf+1):(2*nf)],0))
  nu2_sd = sqrt(pmax(res[(2*nf+1):(3*nf)],0))
  nu3_sd = sqrt(pmax(res[(3*nf+1):(4*nf)],0))
  nu4_sd = sqrt(pmax(res[(4*nf+1):(5*nf)],0))

  if(length(model0)>0) {
      rr = addmom(nu1_sd,model0$nu1_sd,"nu1")
      rr = addmom(nu4_sd,model0$nu2_sd,"nu4",rr)
      rr = addmom(nu2_sd,model0$eps2s_sd,"eps2",rr)
      rr = addmom(nu3_sd,sqrt(1-model0$rho32s^2)*model0$eps3s_sd,"eps3",rr)
      rr = addmom(BEsd,model0$Esd,"Esd",rr)
            
      print(ggplot(rr,aes(x=val2,y=val1,color=type)) + geom_point() + facet_wrap(~name,scale="free") + theme_bw() + geom_abline(linetype=2))

  }

  
  return(list(r1=r1,r4=r4,rt=rt,BEsd=BEsd,nu1_sd=nu1_sd,nu2_sd=nu2_sd,nu4_sd=nu4_sd,nu3_sd=nu3_sd,fit1=obj1,res2=as.numeric(YY - XX %*% res),fitted=as.numeric(XX %*% res),dep=as.numeric(YY)))
}


# extract Var(alpha|l1,l2) and Var(epsion|l)
m4.mini.getvar.stayers.unc2 <- function(sdata,r1,r4,rt,weights=c(),model0=c(),diff=FALSE) {
  
  # USE MOVERS TO GET Var(epsilon|l1)
  setkey(sdata,j1)
  YY1 = sdata[, var(y1)                 ,j1][,V1]   # YY1,YY5,YY6,YY7, YY5,YY2,YY8,YY9, YY6,YY8,YY3,YY0, YY7,YY9,YY0,YY4
  YY2 = sdata[, var(y2)                 ,j1][,V1]
  YY3 = sdata[, var(y3)                 ,j1][,V1]
  YY4 = sdata[, var(y4)                 ,j1][,V1]
  YY5 = sdata[, cov(y1,y2)              ,j1][,V1]
  YY6 = sdata[, cov(y1,y3)              ,j1][,V1]
  YY7 = sdata[, cov(y1,y4)              ,j1][,V1]
  YY8 = sdata[, cov(y2,y3)              ,j1][,V1]
  YY9 = sdata[, cov(y2,y4)              ,j1][,V1]
  YY0 = sdata[, cov(y3,y4)              ,j1][,V1]
  W   = sdata[,.N,j1][,N]
  
  nf  = max(sdata$j1)
  XX1 = array(0,c(nf,5*nf))
  XX2 = array(0,c(nf,5*nf))
  XX3 = array(0,c(nf,5*nf))
  XX4 = array(0,c(nf,5*nf))
  XX5 = array(0,c(nf,5*nf))
  XX6 = array(0,c(nf,5*nf))
  XX7 = array(0,c(nf,5*nf))
  XX8 = array(0,c(nf,5*nf))
  XX9 = array(0,c(nf,5*nf))
  XX0 = array(0,c(nf,5*nf))
  L = nf
  
  for (l1 in 1:nf) {
    ll = l1
    Ls = 0; 
    
    # the Var(alpha|l1,l2)
    XX1[ll,ll] = 1
    XX2[ll,ll] = 1
    XX3[ll,ll] = 1
    XX4[ll,ll] = 1
    XX5[ll,ll] = 1
    XX6[ll,ll] = 1
    XX7[ll,ll] = 1
    XX8[ll,ll] = 1
    XX9[ll,ll] = 1
    XX0[ll,ll] = 1
    Ls= Ls+L
    
    # Var(nu_1|l)
    XX1[ll,Ls+l1] = 1;
    Ls= Ls+L
    
    # Var(nu_2|l)
    XX1[ll,Ls+l1] = r1^2;
    XX2[ll,Ls+l1] = 1;
    XX3[ll,Ls+l1] = rt^2
    XX4[ll,Ls+l1] = r4^2 * rt^2 # here should r4!!!!! @fixme
    XX5[ll,Ls+l1] = r1
    XX6[ll,Ls+l1] = r1*rt
    XX7[ll,Ls+l1] = r1*rt*r4
    XX8[ll,Ls+l1] = rt
    XX9[ll,Ls+l1] = rt*r4
    XX0[ll,Ls+l1] = rt^2*r4 
    Ls= Ls+L
    
    # Var(nu_3|l)
    XX3[ll,Ls+l1] = 1
    XX4[ll,Ls+l1] = r4^2
    XX0[ll,Ls+l1] = r4
    Ls= Ls+L
    
    # Var(nu_4|l)
    XX4[ll,Ls+l1] = 1
  }
  
  XX = rbind(XX1,XX2,XX3,XX4,XX5,XX6,XX7,XX8,XX9,XX0)
  YY =     c(YY1,YY2,YY3,YY4,YY5,YY6,YY7,YY8,YY9,YY0)

  if (length(weights)==0) {
    weights =  c(W,W,W,W,W,W,W,W,W,W)
  } 
  
  if (diff==TRUE) {
    D = t(array(c(-1,1,0,0,
                0,-1,1,0,
                0,0,-1,1),c(4,3)))
    
    DD = kronecker(D,D)
    DD = kronecker(DD,diag(10))
    
    XX = rbind(XX1,XX5,XX6,XX7, XX5,XX2,XX8,XX9, XX6,XX8,XX3,XX0, XX7,XX9,XX0,XX4)
    YY =     c(YY1,YY5,YY6,YY7, YY5,YY2,YY8,YY9, YY6,YY8,YY3,YY0, YY7,YY9,YY0,YY4)

    XX = DD %*% XX
    YY = as.numeric(DD %*% YY)
    XX = XX[,11:50]
    
    weights =  c(W,W,W,W,W,W,W,W,W)
  }
  
  
  fit1 = lm.wfitnn(XX,YY,weights)  
  res  = fit1$solution
  obj1 = t( YY - XX %*% res ) %*% diag(weights) %*% ( YY - XX %*% res )
  
  res2 = res
  if (diff==TRUE) {
    res2 = c(rep(0,10),res)
  }
  
  BEsd   = sqrt(pmax(res2[(0*nf+1):(1*nf)],0))
  nu1_sd = sqrt(pmax(res2[(1*nf+1):(2*nf)],0))
  nu2_sd = sqrt(pmax(res2[(2*nf+1):(3*nf)],0))
  nu3_sd = sqrt(pmax(res2[(3*nf+1):(4*nf)],0))
  nu4_sd = sqrt(pmax(res2[(4*nf+1):(5*nf)],0))
  
  if(length(model0)>0) {
    rr = addmom(nu1_sd,model0$nu1_sd,"nu1")
    rr = addmom(nu4_sd,model0$nu2_sd,"nu4",rr)
    rr = addmom(nu2_sd,model0$eps2s_sd,"eps2",rr)
    rr = addmom(nu3_sd,sqrt(1-model0$rho32s^2)*model0$eps3s_sd,"eps3",rr)
    rr = addmom(BEsd,model0$Esd,"Esd",rr)
    
    print(ggplot(rr,aes(x=val2,y=val1,color=type)) + geom_point() + facet_wrap(~name,scale="free") + theme_bw() + geom_abline(linetype=2))
    
  }
  
  return(list(r1=r1,r4=r4,rt=rt,BEsd=BEsd,nu1_sd=nu1_sd,nu2_sd=nu2_sd,nu4_sd=nu4_sd,nu3_sd=nu3_sd,fit1=obj1,res2=as.numeric(YY - XX %*% res),weights=weights,fitted=as.numeric(XX %*% res),dep=as.numeric(YY)))
}

# extract Var(alpha|l1,l2) and Var(epsion|l)
m4.mini.getvar.movers.unc <- function(jdata,r1,r4) {
  
  # USE MOVERS TO GET Var(epsilon|l1)
  setkey(jdata,j1,j2)
  YY1 = jdata[, var(y1-r1*y2)           ,list(j1,j2)][,V1]
  YY2 = jdata[, cov(y1-r1*y2, y2)       ,list(j1,j2)][,V1]
  YY3 = jdata[, cov(y1-r1*y2, y3)       ,list(j1,j2)][,V1]
  YY4 = jdata[, cov(y1-r1*y2, y4-r4*y3) ,list(j1,j2)][,V1]
  YY5 = jdata[, cov(y2      , y4-r4*y3) ,list(j1,j2)][,V1]
  YY6 = jdata[, cov(y3      , y4-r4*y3) ,list(j1,j2)][,V1]
  YY7 = jdata[, var(y4-r4*y3)           ,list(j1,j2)][,V1]
  W   = jdata[,.N,list(j1,j2)][,N]
  
  nf  = max(jdata$j1)
  XX1 = array(0,c(nf^2,nf^2+2*nf))
  XX2 = array(0,c(nf^2,nf^2+2*nf))
  XX3 = array(0,c(nf^2,nf^2+2*nf))
  XX4 = array(0,c(nf^2,nf^2+2*nf))
  XX5 = array(0,c(nf^2,nf^2+2*nf))
  XX6 = array(0,c(nf^2,nf^2+2*nf))
  XX7 = array(0,c(nf^2,nf^2+2*nf))
  L = nf
  
  for (l1 in 1:nf) for (l2 in 1:nf) {
    ll = l2 + nf*(l1 -1)
    Ls = 0; 
    
    # the Var(alpha|l1,l2)
    XX1[ll,ll] = (1-r1)^2
    XX2[ll,ll] = (1-r1)
    XX3[ll,ll] = (1-r1)
    XX4[ll,ll] = (1-r1)*(1-r4)
    XX5[ll,ll] = (1-r4)
    XX6[ll,ll] = (1-r4)
    XX7[ll,ll] = (1-r4)^2
    Ls= Ls+L^2
    
    # the var(nu|l)
    XX1[ll,Ls+l1] = 1;Ls= Ls+L
    XX7[ll,Ls+l1] = 1
  }
  
  XX = rbind(XX1,XX2,XX3,XX4,XX5,XX6,XX7)
  YY =     c(YY1,YY2,YY3,YY4,YY5,YY6,YY7)
  
  fit1 = lm.wfitnn(XX,YY,c(W,W,W,W,W,W,W))  
  res  = fit1$solution
  obj1 = t( YY - XX %*% res ) %*% diag(c(W,W,W,W,W,W,W)) %*% ( YY - XX %*% res )
  
  EEsd    = t(sqrt(pmax(array((res)[1:nf^2],c(nf,nf)),0)))
  nu1_sd = sqrt(pmax((res)[(nf^2+1):(nf^2+nf)],0))
  nu2_sd = sqrt(pmax((res)[(nf^2+nf+1):(nf^2+2*nf)],0))
  
  return(list(EEsd=EEsd,nu1_sd=nu1_sd,nu2_sd=nu2_sd,fit1=obj1))
}


m4.mini.test <- function() {
  
  # Estimate full model
  model  = m4.mini.new(10,r1=0.4,r4=0.8,eps_sd_factor = 1)
  NNm    = array(50000/(model$nf^2),c(model$nf,model$nf))
  NNs    = array(100000/model$nf,model$nf)
  jdata  = m4.mini.simulate.movers(model,NNm)
  sdata  = m4.mini.simulate.stayers(model,NNs)
  model1 = m4.mini.estimate(jdata,sdata,model$r1,model$r4,model0 = model);
  catf("r1=%f r4=%f r32=%f\n",model$r1,model$r4,model$eps_cor)
  tmp=m4.mini.vardec(model1,verbose=T);
  model$Ns = NNs; model$Nm = NNm
  tmp=m4.mini.vardec(model,verbose=T);
  
  m4.mini.plot(model1)
  
  # Estimate full model with fixed Bs
  model = m4.mini.new(10,fixb=T,r1=0.3,r4=0.6)
  NNm   = array(40000/(model$nf^2),c(model$nf,model$nf))
  NNs  = array(100000/model$nf,model$nf)
  jdata = m4.mini.simulate.movers(model,NNm)
  sdata = m4.mini.simulate.stayers(model,NNs)
  model1=m4.mini.estimate(jdata,sdata,model$r1,model$r4,model0 = model,fixb=T);
  catf("r1=%f r4=%f r32=%f\n",model$r1,model$r4,model$eps_cor)
  
  # Estimate linear model
  model = m4.mini.new(10,linear = T)
  NNm   = array(20000/(model$nf^2),c(model$nf,model$nf))
  NNs  = array(1000000/model$nf,model$nf)
  jdata = m4.mini.simulate.movers(model,NNm)
  sdata = m4.mini.simulate.stayers(model,NNs)
  model1 = m4.minilin.estimate(jdata,sdata,model$r1,model$r4,model0=model)
  #m4.mini.vardec(model1)
  
  # do a grid on r1 and check fit
  dd = expand.grid(r1=seq(0,0.99,l=25),r4=seq(0,0.99,l=25),fit=Inf)
  for (i in 1:nrow(dd)) {
    model1 = m4.minilin.estimate(jdata,sdata,dd$r1[i],dd$r4[i])
    dd$fit1[i] = model1$fit1
    dd$fit2[i] = model1$fit2
    dd$fit3[i] = model1$fit3
    vv = m4.mini.vardec(model1)
    dd$v_psi[i] = vv$v_psi
    dd$v_alpha[i] = vv$v_alpha
    dd$v_2cov[i] = vv$v_2cov
    dd$v_eps[i] = vv$v_eps
    catf("[%i/%i] r1=%f r4=%f fit3=%f \n",i,nrow(dd),dd$r1[i],dd$r4[i],model1$fit3)
  }  
  ggplot(dd,aes(x=rho,y=fit3)) + geom_point() + geom_line() + 
    geom_vline(xintercept=model$r1,color="red",linetype=2) + theme_bw() 

  wireframe(fit3 ~ r1 + r4,dd)
  
  # testing
  
  model1 = m4.mini.getvar(jdata,sdata,model,r1=model$r1,r4=model$r4)
  plot(model$nu1_sd,model1$nu1_sd)
  
  
  m4.mini.getvar.stayers
  
  
  # load the data
  load("../figures/src/mini4p-linear-10.dat",verbose=T)
  m4.mini.vardec(model1)
  m4.mini.vardec(model2)
  
  
  Dmat = diag(c(1059.46526020664,2197.73648768006,5232.6060981711,5986.39143708911,4619.8257387377,2956.56357907282,5872.3510132438,5735.03707424825,2926.84811684315,1496.80395221145))
  Amat = diag(c(1.000000e+00,3.296800e-01,5.575260e-02,1.301464e-01,8.098885e-03,8.432904e-05,3.076671e-03,1.026615e-01,1.911045e-01,2.350272e+00))
  dvec = c(14.012253278123,-5.74691605253042,-63.4752811846924,-41.8857494513356,15.5578043209264,-20.73392047329,-70.9155514168175,-61.0727647947296,8.54506221361423,57.6462542228647)
  bvec = rep(0,10)
  solve.QP(Dmat,as.numeric(dvec),Amat,bvec,meq = 0,factorized = F)  
  
  LowRankQP(Dmat,dvec,array(0,c(1,10)),c(0),uvec=rep(10000,10))
  
  require(LowRankQP)
  
  load("../quadprog_text.dat",verbose=T)
  solve.QP(Dmat,dvec,Amat,bvec)
  qprog(Dmat,dvec,Amat,bvec)  

  # ========= REESTIMATING THE LINEAR MODEL ================ #
  load("../figures/src/mini4p-linear-10.dat",verbose=T)
  NNm   = model1_atm$Nm
  NNs   = model1_atm$Ns
  jdata = m4.mini.simulate.movers(model1_atm,NNm)
  sdata = m4.mini.simulate.stayers(model1_atm,NNs)
  model1 = m4.minilin.estimate(jdata,sdata,model1_atm$r1,model1_atm$r4,model0=model1_atm)
  
  dd = expand.grid(r1=seq(0,0.99,l=25),r4=seq(0,0.99,l=25),fit=Inf)
  for (i in 1:nrow(dd)) {
    model1 = m4.minilin.estimate(jdata,sdata,dd$r1[i],dd$r4[i])
    dd$fit1[i] = model1$fit1
    dd$fit2[i] = model1$fit2
    dd$fit3[i] = model1$fit3
  }
  
  dd$fit1b = pmin(dd$fit1,median(dd$fit1))
  wireframe(-fit1b~r1+r4,dd)
  dd$fit3b = pmin(dd$fit3,median(dd$fit3))
  wireframe(-fit3b~r1+r4,dd)

  g1=ggplot(dd,aes(x=r1,group=r4,y=fit1)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r1,linetype=2,color="red") + geom_vline(xintercept=dd$r1[which.min(dd$fit1)],linetype=2,color="blue")
  g2=ggplot(dd,aes(x=r4,group=r1,y=fit1)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r4,linetype=2,color="Red") + geom_vline(xintercept=dd$r4[which.min(dd$fit1)],linetype=2,color="blue")
  g3=ggplot(dd,aes(x=r1,group=r4,y=fit3b)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r1,linetype=2,color="red") + geom_vline(xintercept=dd$r1[which.min(dd$fit1)],linetype=2,color="blue")
  g4=ggplot(dd,aes(x=r4,group=r1,y=fit3b)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r4,linetype=2,color="Red") + geom_vline(xintercept=dd$r4[which.min(dd$fit1)],linetype=2,color="blue")
  blm:::multiplot(g1,g2,g3,g4,cols=2)
  
  
  dd[which.min(dd$fit1),]

  # ========= REESTIMATING THE FIXB MODEL ================ #
  load("../figures/src/mini4p-fixb-10.dat",verbose=T)
  model0 = model1_atm
  NNm    = 20*model0$Nm
  NNs    = 20*model0$Ns
  jdata  = m4.mini.simulate.movers(model0,NNm)
  sdata  = m4.mini.simulate.stayers(model0,NNs)
  model1 = m4.mini.estimate(jdata,sdata,model0$r1,model0$r4,model0=model0,fixb=T,norm=1,alpha=0)
  
  dd = expand.grid(r1=seq(0,0.99,l=25),r4=seq(0,0.99,l=25),fit=Inf)
  pb <- txtProgressBar(min = 0, max = nrow(dd), style = 3)
  for (i in 1:nrow(dd)) {
    model1 = m4.mini.estimate(jdata,sdata,dd$r1[i],dd$r4[i],fixb=T)
    dd$fit1[i] = model1$fit1
    dd$fit2[i] = model1$fit2
    dd$fit3[i] = model1$fit3
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  dd$fit1b = pmin(dd$fit1,median(dd$fit1))
  wireframe(-fit1b~r1+r4,dd)
  dd$fit3b = pmin(dd$fit3,median(dd$fit3))
  wireframe(-fit3b~r1+r4,dd)

  g1=ggplot(dd,aes(x=r1,group=r4,y=fit1)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r1,linetype=2,color="red") + geom_vline(xintercept=dd$r1[which.min(dd$fit1)],linetype=2,color="blue")
  g2=ggplot(dd,aes(x=r4,group=r1,y=fit1)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r4,linetype=2,color="Red") + geom_vline(xintercept=dd$r4[which.min(dd$fit1)],linetype=2,color="blue")
  g3=ggplot(dd,aes(x=r1,group=r4,y=fit3b)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r1,linetype=2,color="red") + geom_vline(xintercept=dd$r1[which.min(dd$fit3)],linetype=2,color="blue")
  g4=ggplot(dd,aes(x=r4,group=r1,y=fit3b)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r4,linetype=2,color="Red") + geom_vline(xintercept=dd$r4[which.min(dd$fit3)],linetype=2,color="blue")
  blm:::multiplot(g1,g2,g3,g4,cols=2)
    
  # ========= REESTIMATING THE FIXB MODEL ================ #
  load("../figures/src/mini4p-allb-10.dat",verbose=T)
  model0 = model1_atm
  NNm    = 20*model0$Nm
  NNs    = model0$Ns
  jdata  = m4.mini.simulate.movers(model0,NNm)
  sdata  = m4.mini.simulate.stayers(model0,NNs)
  model1 = m4.mini.estimate(jdata,sdata,model0$r1,model0$r4,model0=model0,fixb=F,alpha=1)
  
  dd = expand.grid(r1=seq(0,0.99,l=25),r4=seq(0,0.99,l=25),fit=Inf)
  pb <- txtProgressBar(min = 0, max = nrow(dd), style = 3)
  for (i in 1:nrow(dd)) {
    model1 = m4.mini.estimate(jdata,sdata,dd$r1[i],dd$r4[i],fixb=F,alpha=1)
    dd$fit1[i] = model1$fit1
    dd$fit2[i] = model1$fit2
    dd$fit3[i] = model1$fit3
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  dd$fit1b = pmin(dd$fit1,median(dd$fit1))
  wireframe(-fit1b~r1+r4,dd)
  dd$fit3b = pmin(dd$fit3,median(dd$fit3))
  wireframe(-fit3b~r1+r4,dd)
  
  g1=ggplot(dd,aes(x=r1,group=r4,y=fit1)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r1,linetype=2,color="red") + geom_vline(xintercept=dd$r1[which.min(dd$fit1)],linetype=2,color="blue")
  g2=ggplot(dd,aes(x=r4,group=r1,y=fit1)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r4,linetype=2,color="red") + geom_vline(xintercept=dd$r4[which.min(dd$fit1)],linetype=2,color="blue")
  g3=ggplot(dd,aes(x=r1,group=r4,y=fit3b)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r1,linetype=2,color="red") + geom_vline(xintercept=dd$r1[which.min(dd$fit3)],linetype=2,color="blue")
  g4=ggplot(dd,aes(x=r4,group=r1,y=fit3b)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r4,linetype=2,color="red") + geom_vline(xintercept=dd$r4[which.min(dd$fit3)],linetype=2,color="blue")
  blm:::multiplot(g1,g2,g3,g4,cols=2)
  
  dd[which.min(dd$fit),]
    
  # -------------- STAYERS --------------------- #
  dd = expand.grid(r1=seq(0,0.99,l=25),r4=seq(0,0.99,l=25),fit=Inf)
  pb <- txtProgressBar(min = 0, max = nrow(dd), style = 3)
  dd[,paste("res",1:7,sep='')]=0
  for (i in 1:nrow(dd)) {
    model1 = m4.mini.getvar.stayers.unc(sdata,dd$r1[i],dd$r4[i])
    dd$fit[i] = model1$fit1
    dd[i,paste("res",1:7,sep='')]=model1$res
    setTxtProgressBar(pb, i)
  }
  close(pb)
  g1=ggplot(dd,aes(x=r1,group=r4,y=fit)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r1,linetype=2,color="red") + geom_vline(xintercept=dd$r1[which.min(dd$fit)],linetype=2,color="blue")
  g2=ggplot(dd,aes(x=r4,group=r1,y=fit)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r4,linetype=2,color="red") + geom_vline(xintercept=dd$r4[which.min(dd$fit)],linetype=2,color="blue")
  blm:::multiplot(g1,g2,cols=2)
  dd[which.min(dd$fit),]
  
  # -------------- MOVERS --------------------- #
  dd = expand.grid(r1=seq(0,0.99,l=25),r4=seq(0,0.99,l=25),fit=Inf)
  pb <- txtProgressBar(min = 0, max = nrow(dd), style = 3)
  for (i in 1:nrow(dd)) {
    model1 = m4.mini.getvar.movers.unc(jdata,dd$r1[i],dd$r4[i])
    dd$fit[i] = model1$fit1
    setTxtProgressBar(pb, i)
  }
  close(pb)
  g1=ggplot(dd,aes(x=r1,group=r4,y=fit)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r1,linetype=2,color="red") + geom_vline(xintercept=dd$r1[which.min(dd$fit)],linetype=2,color="blue")
  g2=ggplot(dd,aes(x=r4,group=r1,y=fit)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r4,linetype=2,color="red") + geom_vline(xintercept=dd$r4[which.min(dd$fit)],linetype=2,color="blue")
  blm:::multiplot(g1,g2,cols=2)
  
  # ------- SIMULATE FROM MIXTURE MODEL AND TRY TO EXTRACT rhos ----------- #
  load("../figures/src/em-endo-levelmodel-6x10.dat",verbose=T)
  model0 = res_het
}

# use ks=c(1,1) for normal and c(0.5,0.1) for high kurtosis
rnorms <- function(n,ks) {
  s2= ks[1]
  lambda = ks[2]
  s1 = sqrt(1- (1-lambda)*s2^2)/sqrt(lambda)
  s1 = sqrt(1- (1-lambda)*s2^2)/sqrt(lambda)
  co = runif(n)<lambda
  X  = rnorm(n)*( co*s1 + (1-co)*s2)
  return(X)
}



test.kurtosis <- function() {
  # simulate from a function with given kurosis and skewness and mean 0 and variance 1
  
  # 2 point normal mixture lamba, m1, m2, s1, s2
  # but mean 0 and variance 1 so we have constraints!
  # take (m1,lambda,s1,s2)
  
  sim.mix <- function(lambda,s2,n=1000) {
    s1 = sqrt(1- (1-lambda)*s2^2)/sqrt(lambda)
    co = runif(n)<lambda
    X  = rnorm(n)*( co*s1 + (1-co)*s2)
    return(c(var(X),skewness(X),kurtosis(X)))
  }
  
  rnorms <- function(n,lambda,s2) {
    s1 = sqrt(1- (1-lambda)*s2^2)/sqrt(lambda)
    s1 = sqrt(1- (1-lambda)*s2^2)/sqrt(lambda)
    co = runif(n)<lambda
    X  = rnorm(n)*( co*s1 + (1-co)*s2)
    return(X)
  }
  
  ff  <- function(theta) sum( (sim.mix(theta[1],theta[2],theta[3],theta[4]) - c(1,2,10))^2)
  ff2 <- function(theta) sim.mix(theta[1],theta[2],theta[3],theta[4])
  res = optim(c(0,0.1,1,10),ff)
  ff2(res$par)
  
  
}



test.getrhos.unc2 <- function() {
  load("../figures/src/mini4p-linear-10.dat",verbose=T)
  model0 = model1_atm
  NNm   = model0$Nm
  NNs   = 15*model0$Ns
  jdata = m4.mini.simulate.movers(model0,NNm)
  
  model0$r1=0.6
  model0$r4=0.7
  model0$rho32s=0.4
  sdata = m4.mini.simulate.stayers(model0,NNs)

  # run at true rhos 
  res = m4.mini.getvar.stayers.unc2(sdata,model0$r1,model0$r4,model0$rho32s,model0=model0)

  nr = 15
  dd = expand.grid(r1=seq(0,0.99,l=nr),r4=seq(0,0.99,l=nr),rt=seq(0,0.99,l=nr),fit=Inf)
  pb <- txtProgressBar(min = 0, max = nrow(dd), style = 3)
  dd[,paste("res",1:10,sep='')]=0
  for (i in 1:nrow(dd)) {
    model1 = m4.mini.getvar.stayers.unc2(sdata,dd$r1[i],dd$r4[i],dd$rt[i])
    dd$fit[i] = model1$fit1
    dd[i,paste("res",1:10,sep='')]=model1$res
    setTxtProgressBar(pb, i)
  }
  close(pb)
  g1=ggplot(dd,aes(x=r1,y=fit)) + geom_point() + theme_bw() + geom_vline(xintercept=model0$r1,linetype=2,color="red") + geom_vline(xintercept=dd$r1[which.min(dd$fit)],linetype=2,color="blue")
  g2=ggplot(dd,aes(x=r4,y=fit)) + geom_point() + theme_bw() + geom_vline(xintercept=model0$r4,linetype=2,color="red") + geom_vline(xintercept=dd$r4[which.min(dd$fit)],linetype=2,color="blue")
  g3=ggplot(dd,aes(x=rt,y=fit)) + geom_point() + theme_bw() + geom_vline(xintercept=model0$rho32s,linetype=2,color="red") + geom_vline(xintercept=dd$rt[which.min(dd$fit)],linetype=2,color="blue")
  multiplot(g1,g2,g3,cols=3)
  dd[which.min(dd$fit),]
  
  # testing strong kurtosis effects
  load("../figures/src/mini4p-fixb-10.dat",verbose=T)
  model0 = model1_atm
  NNm   = model0$Nm
  NNs   = 5*model0$Ns
  jdata = m4.mini.simulate.movers(model0,NNm)
  
  model0$r1=0.4
  model0$r4=0.6
  model0$rho32s=0.3
  sdata = m4.mini.simulate.stayers(model0,NNs,ks = c(0.5,0.1))
  sdata = m4.mini.simulate.stayers(model0,NNs,ks = c(1,1))
  
  # estimate mixture of normal
  modelr = em.endo.level.new(6,10)
  ctrl   = em.control(nplot=10,check_lik=F,fixb=T,est_rho=T)
  res    = em.endo.level.rhoext.stayers.unc(jdata,sdata,modelr,ctrl)
  

}

test.simulate.fromest <- function() {
  
  source("R/utils.R")
  path_archive = "../figures/res/"
  arch_dyn = "res-2003-dynamic.dat"
  
  archive.list(arch_dyn)
  mini_model       = archive.get("mini_model",arch_dyn)
  cstats             = archive.get("cstats",arch_dyn)
  
  sdata = m4.mini.simulate.stayers(mini_model,mini_model$Ns)
  sdata[, j2:=j1]
  sdata = m4.mini.impute.stayers(mini_model,sdata)

  jdata = m4.mini.simulate.movers(mini_model,mini_model$Nm)
  jdata = m4.mini.impute.movers(mini_model,jdata)
  
  adata = rbind(
            sdata[,list(y1=y1_imp,y2=y2_imp,y3=y3_imp,y4=y4_imp,j1,j2,k_imp,move=FALSE)],
            jdata[,list(y1=y1_imp,y2=y2_imp,y3=y3_imp,y4=y4_imp,j1,j2,k_imp,move=TRUE)])
  
  save(adata,file="~/Dropbox/sortingwithlda/data/simdata-dyn.dat")
  
  load("~/Dropbox/sortingwithlda/data/simdata-dyn.dat")
  
  # check the ordering of the clusters
  summary(lm(y2~ 0 + factor(j1),adata))
  
  # check the effect of move
  summary(lm(y2~ move + k_imp + factor(j1),adata))
  summary(lm(y2~ move + k_imp:factor(j1),adata))
  
  # check the effect of the firm in period 2
  summary(lm(y2~ k_imp + factor(j1) + factor(j2),adata[move==TRUE]))
  summary(lm(y2~ k_imp:factor(j1) + factor(j2),adata[move==TRUE]))
  
  # check the effect on y3
  summary(lm(y3~ k_imp:factor(j2) + factor(j1),adata[move==TRUE]))
  summary(lm(y3~ k_imp:factor(j2) + factor(j1) + y2,adata[move==TRUE]))

}

test.simulate.from_mc_est <- function() {
  
  source("R/utils.R")
  path_archive = "../figures/res/"
  arch_dyn = "res-2003-dynamic.dat"
  
  archive.list(arch_dyn)
  mini_models       = archive.get("mini_bs",arch_dyn)
  cstats            = archive.get("cstats",arch_dyn)
  
  getc <- function(fit,i) {
    rr1 = data.frame(coefficients(summary(fit)))
    rr1$name = paste(formula(fit)[2:3],collapse = " ~ ")
    rr1$i=i
    rr1$variable = rownames(rr1)
    rownames(rr1)<-NULL
    return(rr1)
  }    

  rrr = data.frame()
  for (i in 1:100) {
    mini_model = mini_models[[i]]
    sdata = m4.mini.simulate.stayers(mini_model,mini_model$Ns)
    sdata[, j2:=j1]
    sdata = m4.mini.impute.stayers(mini_model,sdata)
    
    jdata = m4.mini.simulate.movers(mini_model,mini_model$Nm)
    jdata = m4.mini.impute.movers(mini_model,jdata)
    
    adata = rbind(
      sdata[,list(y1=y1_imp,y2=y2_imp,y3=y3_imp,y4=y4_imp,j1,j2,k_imp,move=FALSE)],
      jdata[,list(y1=y1_imp,y2=y2_imp,y3=y3_imp,y4=y4_imp,j1,j2,k_imp,move=TRUE)])
    
    fit = lm(y2~ 0 + factor(j1),adata)
    rrr = rbind(rrr,getc(fit,i))
    fit = lm(y2~ move + k_imp + factor(j1),adata)
    rrr = rbind(rrr,getc(fit,i))
    fit = lm(y2~ move + k_imp:factor(j1),adata)
    rrr = rbind(rrr,getc(fit,i))
    fit = lm(y2~ k_imp + factor(j1) + factor(j2),adata[move==TRUE])
    rrr = rbind(rrr,getc(fit,i))
    fit = lm(y2~ k_imp:factor(j1) + factor(j2),adata[move==TRUE])
    rrr = rbind(rrr,getc(fit,i))
    fit = lm(y3~ k_imp:factor(j2) + factor(j1),adata[move==TRUE])
    rrr = rbind(rrr,getc(fit,i))
    fit = lm(y3~ k_imp:factor(j2) + factor(j1) + y2,adata[move==TRUE])
    rrr = rbind(rrr,getc(fit,i))
    fit = lm(y3~ k_imp+ factor(j2) + factor(j1) + y2,adata[move==TRUE])
    rrr = rbind(rrr,getc(fit,i))
    catf("done with %i\n",i)
  }
  
  rrr  = data.table(rrr)
  rrr2 = rrr[, list(m=mean(Estimate),q1=quantile(Estimate,0.025),q2=quantile(Estimate,0.975),sd=sd(Estimate)),list(name,variable)] 
  
  # plotting bootstrapped parameters
  
  
  rr= ldply(mini_models,function(x) data.frame(k=1:10,B=x$B1))
  ggplot(rr,aes(x=B)) + geom_histogram() + facet_wrap(~k,scale="free")
  rr= ldply(mini_models,function(x) { 
    r = expand.grid(k1=1:10,k2=1:10); r$value=c(x$EEsd);r}
    )
  rr = data.table(rr)
  ggplot(rr[k1%in%c(2,4,6,8)],aes(x=value)) + geom_histogram() + facet_grid(k1~k2)
  
  # append replication
    
}
