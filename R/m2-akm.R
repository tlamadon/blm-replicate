#' file to run AKM

akm.extract.getjdata <- function(sim.data,year1=2002,year2=2004) {
  jdata.sim = sim.data[ (t==year1) | (t==year2)]
  setkey(jdata.sim,wid,t)
  jdata.sim = jdata.sim[ , list(j1 = l[1],j2=l[2],y1=y[1],y2=y[2],
                                f1=fid[1],f2=fid[2],x=x[1]),wid]
  return(jdata.sim)
}

#' Function to create CSR sparse matrices
#' @export
sparseMatrix2 <- function(i,j,v,dim) {
  
  A = matrix(1,2,2)
  A = as.matrix.csr(A)
  
  A@dimension = as.integer(dim)
  
  # get the indices
  dd = data.table(i=i,j=j,v=v)
  setkey(dd,i,j)
  dd[,c:=1:.N]
  
  ia = dd[,.N,i]
  # add empty rows
  ia2 = rep(0,dim[1])
  ia2[ia$i]=ia$N
  ia = cumsum(ia2)
  ia =c(ia,ia[length(ia)]+1)
  
  A@ra = dd$v
  A@ja = as.integer(dd$j)
  A@ia = as.integer(ia)
  
  return(A)
}

#' Compute the firm fix effects using the movers. This first 
#' extracts the largest connected set among the firms, then uses
#' pairs in a sparse regression.
#' 
#' @export
m2.fe.firms <-function(jdata,fix_sd=NA) {
  
  # start by extracting the largest connected set
  jdata = jdata[,j1:=f1][,j2:=f2]
  f1s   = get.largest.conset(jdata)
  if (is.numeric(jdata$j1)) f1s = as.integer(f1s);
  
  jdata = jdata[f1  %in% f1s][f2 %in% f1s]
  
  # reindex the firms (@e: why?)
  fids = data.table(f1=f1s,nfid=1:length(f1s))
  setkey(fids,f1)
  setkey(jdata,f1)
  jdata[,j1 := fids[jdata,nfid]]
  setkey(jdata,f2)
  jdata[,j2 := fids[jdata,nfid]]
  
  nf = pmax(max(jdata$j1),max(jdata$j2))
  # _________   step 1, get the psi ___________
  #compute the means (across individuals who make the movement from firm j1 to firm j2 (mostly 1))
  #dd = jdata[,list(m1 = mean(y2 - y1),.N),list(j1,j2)]
  #N = nrow(dd)
  #XX1 = matrix(0,N,nf)
  #XX1[1:N + N*(dd$j1-1)]=1
  #XX2 = matrix(0,N,nf)
  #XX2[1:N + N*(dd$j2-1)]=1
  #XX  = as.matrix.csr(XX2) - as.matrix.csr(XX1)
  #fit = SparseM::slm.fit(XX[,2:nf],dd$m1,nnzlmax=1e8,tmpmax=1e8,nsubmax=1e8)$coefficients
  
  dd = jdata[,list(m1 = mean(y2 - y1),.N),list(j1,j2)]
  N = nrow(dd)
  dd[,v1:=-1*(j1>1)][,v2:=1*(j2>1)]
  dd[,c1:=ifelse(j1==1,1,j1-1)][,c2:=ifelse(j2==1,1,j2-1)]
  XX = sparseMatrix2(1:N,dd$c1,dd$v1,c(N,nf-1)) + sparseMatrix2(1:N,dd$c2,dd$v2,c(N,nf-1))
  fit = SparseM::slm.fit(XX,dd$m1,nnzlmax=1e8,tmpmax=1e8,nsubmax=1e8)$coefficients
  
  Psi  = as.numeric(c(0,c(fit)))
  fids[,psi := Psi[nfid]]
  
  # attach the number of movers for each firm
  tmp = rbind(jdata[,.N,list(f=f1)],jdata[,.N,list(f=f2)])
  tmp = tmp[,sum(N),f]
  setkey(tmp,f)
  setkey(fids,f1)
  fids[, nmovers:= tmp[fids,V1]]

  # ________ get variances __________
  V_eps1 = jdata[,var(y2 -y1 - Psi[j2]+Psi[j1])]/2
  
  return(list(var_eps=V_eps1,fe=fids))
}


#' Compute the fixed-effect estimator
#' of the wokrer and firm decomposition
#' 
#' @export
m2.fe.all <-function(sim) {
  
  # we start by extracting the firm FE
  # on the connected set
  data.fe = m2.fe.firms(sim$jdata)$fids
  
  ad = rbind(sim$sdata[,list(f1,y1,f2,y2)],sim$jdata[,list(f1,y1,f2,y2)])
  ad = ad[ f1 %in% data.fe$f1 ][ f2 %in% data.fe$f2 ]
  
  # then we merge back the firm FE into period 1
  setkey(data.fe,f1)
  setkey(ad,f1)
  ad[,psi1 := data.fe[ad]]
  setkey(ad,f2)
  ad[,psi2 := data.fe[ad]]

  # then compute worker FE
  ad[,alpha := 0.5*(y1 - psi + y2-psi2)]
  
  # finally we compute the decomposition
  cc = ad[,cov(cbind(y1,x,psi1,y1-psi1-x))]  
  
  get.stats <- function(cc) {
    r=list()
    den = cc[2,2] + cc[3,3] + 2 * cc[2,3]    
    r$cor_kl = round(cc[2,3]/sqrt(cc[2,2]*cc[3,3]),4)
    r$cov_kl = 2*round(cc[2,3]/den,4)
    r$var_k  = round(cc[2,2]/den,4)
    r$var_l  = round(cc[3,3]/den,4)
    r$rsq    = round((cc[1,1] - cc[4,4])/cc[1,1],4)
    return(r)
  }
  
  stats = get.stats(cc)
  return(list(cc=cc,stats=stats))
}


#' Extracts the largest connected set from data on 
#' movers
#' @export
get.largest.conset <- function(jdata) {
  # combine th 2 directions
  jdata2 = rBind(jdata[,list(j1,j2)],jdata[,list(j1=j2,j2=j1)])
  AD = pmin(acast(jdata2[,.N,list(j1,j2)],value.var = "N",j1~j2,drop=FALSE,fill=0),1)
  # compute connected sets
  cs = conComp(AD,2)
  # extract largest
  cs = data.table(f1=names(cs),set=as.numeric(cs))
  cs.size = cs[,.N,set][order(-N)]
  return(cs[  set == cs.size$set[1]  , f1])
}



