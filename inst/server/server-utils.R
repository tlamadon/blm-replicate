# some utility functions on the server

res.load <- function(name) {
  destfile = sprintf("%s/%s.dat",local_opts$wdir,name)
  flog.info("loading %s",destfile)

  # check that dependent resource exists
  if (!file.exists(destfile)) {
    flog.error("Resource %s does not exists, you need to run the code that creates it.",name)
  }

  load(destfile)
  return(value)
}

res.save <- function(name,value) {
  destfile = sprintf("%s/%s.dat",local_opts$wdir,name)
  flog.info("saving to %s",destfile)
  save(value,file=destfile)
}



#' combines the different legs of the bootstraps
m4.getboostrap <- function() {
  res_bs         = rkiv0.load("m4-mixt-d2003-bootstrap") # we do it across bootstraps
  res_bs2        = rkiv0.load("m4-mixt-d2003-bootstrap-leg2") # we do it across bootstraps
  names(res_bs2$mixt_all) = paste(101:(100+length(res_bs2$mixt_all)))
  res_bs2$vdecs$rep = res_bs2$vdecs$rep+100
  res_bs$mixt_all = c(res_bs$mixt_all,res_bs2$mixt_all)
  res_bs$vdecs = data.table(rbind(res_bs$vdecs,res_bs2$vdecs))
  return(res_bs)
}

#' combines the different legs of the bootstraps
m2.getboostrap <- function(lim=200) {
  res_bs       = rkiv0.load("m2-mixt-d2003-bootstrap") # we do it across bootstraps
  res_bs2      = rkiv0.load("m2-mixt-d2003-bootstrap-leg2") # we do it across bootstraps
  names(res_bs2$mixt_all) = paste(101:(100+length(res_bs2$mixt_all)))
  res_bs$mixt_all = c(res_bs$mixt_all,res_bs2$mixt_all)
  res_bs2$vdecs$rep = res_bs2$vdecs$rep + 100
  res_bs$vdecs = data.table(rbind(res_bs$vdecs,res_bs2$vdecs))
  res_bs$vdecs = res_bs$vdecs[rep<lim]
  res_bs$mixt_all = res_bs$mixt_all[paste(1:lim)]
  return(res_bs)
}

get.stats.clusters <- function(data,movers=FALSE,ydep="y1") {
  rr = list()
  rr$nwid       = data[,length(unique(wid))]   # unique worker
  rr$nfid       = data[,length(unique(f1))]

  # creating firm info from 2002
  fdata = data[,list(.N,ind=ind1[1],size1=size1[1],va1=va1[1],Nm=sum(move==TRUE)),f1]

  rr$nfirm_actualsize_ge10       = fdata[N>=10,.N]
  rr$nfirm_actualsize_ge50       = fdata[N>=50,.N]
  rr$nfirm_reportedsize_ge10       = fdata[size1>=10,.N]
  rr$nfirm_reportedsize_ge50       = fdata[size1>=50,.N]
  rr$nfirm_movers_ge1       = fdata[Nm>=1,.N]
  rr$nfirm_movers_ge5       = fdata[Nm>=5,.N]
  rr$nfirm_movers_ge10      = fdata[Nm>=10,.N]
  rr$firm_reportedsize_mean   = fdata[,mean(size1)]
  rr$firm_reportedsize_median = fdata[,median(size1)*1.0]
  rr$firm_actualsize_mean   = fdata[,mean(N)]
  rr$firm_actualsize_median = fdata[,median(N)*1.0]

  rr$firm_reportedsize_median_worker = data[,median(size1)*1.0]
  rr$firm_actualsize_median_worker = data[,median(asize)*1.0]
  rr$firm_mean_log_va = fdata[va1>0,mean(log(va1))]
  rr$firm_var_log_va = fdata[va1>0,var(log(va1))]
  rr$firm_neg_va = fdata[va1<=0,.N]

  rr$worker_share_educ1   = data[,mean(educ==1)]
  rr$worker_share_educ2   = data[,mean(educ==2)]
  rr$worker_share_educ3   = data[,mean(educ==3)]
  rr$worker_var_log_wage  = data[,var(get(ydep))]
  rr$worker_mean_log_wage = data[,mean (get(ydep))]
  rr$worker_skewness_log_wage = data[,skewness(get(ydep))]
  rr$worker_kurtosis_log_wage = data[,kurtosis(get(ydep))]

  rr$worker_share_ind_manu   = data[,mean(ind1=="Manufacturing")]
  rr$worker_share_ind_serv   = data[,mean(ind1=="Services")]
  rr$worker_share_ind_retail   = data[,mean(ind1=="Retail trade")]
  rr$worker_share_ind_cons   = data[,mean(ind1=="Construction etc.")]

  # between firm variance
  rr$between_firm_wage_var = data[, var(fmw)]

  rr$worker_share_age_0_30   = data[,mean(age<=30)]
  rr$worker_share_age_31_50  = data[,mean(age>=31 & age<=50)]
  rr$worker_share_age_51_inf = data[,mean(age>=51)]

  return(rr)
}

compute.mob.matrix <- function(model,nsim=1e7) {
  # we simulate serveral time, and collect based on the realization of epsilon in period 2
  NNm = model$NNm
  NNs = model$NNs
  NNm[!is.finite(NNm)]=0
  NNs[!is.finite(NNs)]=0
  stayer_share = sum(10*NNs)/(sum(10*NNs)+sum(NNm))
  NNs = round(NNs*nsim*stayer_share/sum(NNs))
  NNm = round(NNm*nsim*(1-stayer_share)/sum(NNm))

  # we simulate from the model both movers and stayers
  sdata.sim = m4.mixt.simulate.stayers(model,NNs)
  jdata.sim = m4.mixt.simulate.movers(model,NNm)
  data.sim = rbind(sdata.sim[,list(k,j1,j2=j1,y2,m=0)],jdata.sim[,list(k,j1,j2,y2,m=1)])

  # make sure the cluster are ordered
  #If = rank(data.sim[,mean(y2),j1][order(j1),V1])
  #data.sim[,j1 := If[j1]][,j2:=If[j2]]

  rm(sdata.sim,jdata.sim)

  # compute the residuals
  data.sim[,er:=rank(y2)/.N,list(j1,k)]
  # aggregate clusters
  ff = function(x) {
    if (x<=3) return("k=1..3");
    if (x<=7) return("k=4..7");
    return("k=8..10");
  }
  data.sim[,j1g := ff(j1),j1]
  data.sim[,j2g := ff(j2),j2]

  # code j2=0 when m=0
  data.sim[m==0,j2g:="0"]

  M1 = acast(data.sim[(er<0.1),.N,list(j1g,j2g,m)],j1g~j2g,fill=0,value.var = "N")
  M1 = M1/spread(rowSums(M1),2,4); M1[,1]=1-M1[,1]
  M2 = acast(data.sim[,.N,list(j1g,j2g,m)],j1g~j2g,fill=0,value.var = "N")
  M2 = M2/spread(rowSums(M2),2,4); M2[,1]=1-M2[,1]
  M3 = acast(data.sim[(er>.9),.N,list(j1g,j2g,m)],j1g~j2g,fill=0,value.var = "N")
  M3 = M3/spread(rowSums(M3),2,4); M3[,1]=1-M3[,1]

  d1 = melt(M1,c('j1','j2'))
  d1$cond=1
  d2 = melt(M2,c('j1','j2'))
  d2$cond=2
  d3 = melt(M3,c('j1','j2'))
  d3$cond=3
  return(data.table(rbind(d1,d2,d3)))
}


analysis.dynamic.dec <- function(model) {

  nsim=1e7
  NNm = model$NNm
  NNs = model$NNs
  NNm[!is.finite(NNm)]=0
  NNs[!is.finite(NNs)]=0
  stayer_share = sum(10*NNs)/(sum(10*NNs)+sum(NNm))
  NNs = round(NNs*nsim*stayer_share/sum(NNs))
  NNm = round(NNm*nsim*(1-stayer_share)/sum(NNm))

  # we only use movers here
  jdata.sim = m4.mixt.simulate.movers(model,20*NNm)

  # depdendence effect
  res = list()

  # Y2
  jdata.sim[,v1  := var(y2),list(j1,j2,k)]
  res$pe_y2_c2  = jdata.sim[,mean(v1,na.rm=T)]
  jdata.sim[,m2 := mean(y2),list(j1,j2,k)]
  jdata.sim[,v2 := var(m2),list(j1,k)]
  res$pe_y2_c1  = jdata.sim[,mean(v2,na.rm=T)]
  jdata.sim[,v3 := var(y2), list(j1,k)]
  res$pe_y2_tot = jdata.sim[,mean(v3,na.rm=T)]

  # Y1
  jdata.sim[,v1  := var(y1),list(j1,j2,k)]
  res$pe_y1_c2  = jdata.sim[,mean(v1,na.rm=T)]
  jdata.sim[,m2 := mean(y1),list(j1,j2,k)]
  jdata.sim[,v2 := var(m2),list(j1,k)]
  res$pe_y1_c1  = jdata.sim[,mean(v2,na.rm=T)]
  jdata.sim[,v3 := var(y1), list(j1,k)]
  res$pe_y1_tot = jdata.sim[,mean(v3,na.rm=T)]

  # Y3
  jdata.sim[,v1  := var(y3),list(j1,j2,k)]
  res$pe_y3_c2  = jdata.sim[,mean(v1,na.rm=T)]
  jdata.sim[,m2 := mean(y3),list(j1,j2,k)]
  jdata.sim[,v2 := var(m2),list(j2,k)]
  res$pe_y3_c1  = jdata.sim[,mean(v2,na.rm=T)]
  jdata.sim[,v3 := var(y3), list(j2,k)]
  res$pe_y3_tot = jdata.sim[,mean(v3,na.rm=T)]

  # also compute y1,y4
  jdata.sim[,v1  := var(y4),list(j1,j2,k)]
  res$pe_y4_c2  = jdata.sim[,mean(v1,na.rm=T)]
  jdata.sim[,m2 := mean(y4),list(j1,j2,k)]
  jdata.sim[,v2 := var(m2),list(j2,k)]
  res$pe_y4_c1  = jdata.sim[,mean(v2,na.rm=T)]
  jdata.sim[,v3 := var(y4), list(j2,k)]
  res$pe_y4_tot = jdata.sim[,mean(v3,na.rm=T)]

  # ---------  firm effect -------------
  jdata.sim[,m3 := mean(y2),list(j1,k)]
  jdata.sim[,v3 := var(m3),list(k)]
  res$fe_y2 = jdata.sim[,mean(v3)]
  jdata.sim[,m3 := mean(y3),list(j2,k)]
  jdata.sim[,v3 := var(m3),list(k)]
  res$fe_y3 = jdata.sim[,mean(v3)]
  jdata.sim[,m3 := mean(y1),list(j2,k)]
  jdata.sim[,v3 := var(m3),list(k)]
  res$fe_y1 = jdata.sim[,mean(v3)]
  jdata.sim[,m3 := mean(y4),list(j2,k)]
  jdata.sim[,v3 := var(m3),list(k)]
  res$fe_y4 = jdata.sim[,mean(v3)]

  # ---------  network effect -------------

  # we need to compute the following 2 terms
  # sum_{k'}p(k'|k)E(Y3|k') + sum_{k'}p(k'|k)(E(Y3|k',k)-E(Y3|k'))

  # first term
  jdata.sim[, m1 := mean(y3), list(k,j2)]
  jdata.sim[, a1 := mean(m1), list(k,j1)]
  # last term
  jdata.sim[, m2 := mean(y3), list(k,j1,j2)]
  jdata.sim[, a2 := mean(m2-m1), list(k,j1)]
  # total
  jdata.sim[, m3 := mean(y3), list(k,j1)]
  rr3 = jdata.sim[, { vv = cov(cbind(a1,a2)); data.frame(vot=var(m3),v1=vv[1,1],v2=vv[2,2],cc=vv[1,2],N=.N)},k]
  rr3 = rr3[,list(v1=wt.mean(v1,N),v2=wt.mean(v2,N),cc=wt.mean(cc,N))]
  res$ne_y3_v1  = rr3$v1
  res$ne_y3_v2  = rr3$v2
  res$ne_y3_cov = rr3$cc


  # first term
  jdata.sim[, m1 := mean(y4), list(k,j2)]
  jdata.sim[, a1 := mean(m1), list(k,j1)]
  # last term
  jdata.sim[, m2 := mean(y4), list(k,j1,j2)]
  jdata.sim[, a2 := mean(m2-m1), list(k,j1)]
  # total
  jdata.sim[, m3 := mean(y4), list(k,j1)]
  rr4 = jdata.sim[, { vv = cov(cbind(a1,a2)); data.frame(vot=var(m3),v1=vv[1,1],v2=vv[2,2],cc=vv[1,2],N=.N)},k]
  rr4 = rr4[,list(v1=wt.mean(v1,N),v2=wt.mean(v2,N),cc=wt.mean(cc,N))]
  res$ne_y4_v1  = rr4$v1
  res$ne_y4_v2  = rr4$v2
  res$ne_y4_cov = rr4$cc

  return(data.frame(res))
}

analysis.dynamic.dec.bis <- function(model,nsim=1e7) {

  NNm = model$NNm
  NNs = model$NNs
  NNm[!is.finite(NNm)]=0
  NNs[!is.finite(NNs)]=0
  stayer_share = sum(10*NNs)/(sum(10*NNs)+sum(NNm))
  NNs = round(NNs*nsim*stayer_share/sum(NNs))
  NNm = round(NNm*nsim*(1-stayer_share)/sum(NNm))

  NNm = round(nsim * NNm/sum(NNm))

  # we only use movers here
  jdata.sim = m4.mixt.simulate.movers(model,NNm)
  #sdata.sim = m4.mixt.simulate.stayers(model,NNs)

  # reorder
  #If = rank(sdata.sim[,mean(y2),j1][order(j1),V1])
  #jdata.sim[,j1 := If[j1]][,j2:=If[j2]]

  # ----- Var(Y1), Y2, Y3, Y4 within alpha ------- #
  res = list()
  v_y1_unc = jdata.sim[,v_y1 := var(y1),k][,mean(v_y1)]
  v_y2_unc = jdata.sim[,v_y2 := var(y2),k][,mean(v_y2)]
  v_y3_unc = jdata.sim[,v_y3 := var(y3),k][,mean(v_y3)]
  v_y4_unc = jdata.sim[,v_y4 := var(y4),k][,mean(v_y4)]

  res$v_y1_unc=v_y1_unc
  res$v_y2_unc=v_y2_unc
  res$v_y3_unc=v_y3_unc
  res$v_y4_unc=v_y4_unc
  #jdata.sim[,v_y2 := var(y2),k]
  #jdata.sim[,v_y3 := var(y3),k]
  #jdata.sim[,v_y4 := var(y4),k]

  # -----  Wa Bk of Y1, Y2    -----  #
  jdata.sim[, m1 := mean(y1), list(k,j1)]
  jdata.sim[, v1 := var(m1), list(k)]
  res$y1_Wa_Bk1 = jdata.sim[, mean(v1/v_y1_unc)]

  jdata.sim[, m1 := mean(y2), list(k,j1)]
  jdata.sim[, v1 := var(m1), list(k)]
  res$y2_Wa_Bk1 = jdata.sim[, mean(v1/v_y2_unc)]

  # -------  Wa,k Bk' of Y1,Y2 ------- #
  jdata.sim[, m1 := mean(y1), list(k,j1,j2)]
  jdata.sim[, v1 := var(m1), list(k,j1)]
  res$y1_Wak1_Bk2 = jdata.sim[, mean(v1/v_y1_unc)]

  jdata.sim[, m1 := mean(y2), list(k,j1,j2)]
  jdata.sim[, v1 := var(m1), list(k,j1)]
  res$y2_Wak1_Bk2 = jdata.sim[, mean(v1/v_y2_unc)]

  # ---------- Wa,k,k' of Y1,Y2 -------- #
  jdata.sim[, v1 := var(y1), list(k,j1,j2)]
  jdata.sim[is.na(v1),v1:=0]
  res$y1_Wak1k2 = jdata.sim[, mean(v1/v_y1_unc)]

  jdata.sim[, v1 := var(y2), list(k,j1,j2)]
  jdata.sim[is.na(v1),v1:=0]
  res$y2_Wak1k2 = jdata.sim[, mean(v1/v_y2_unc)]

  # -----  Wa Bk' of Y3, Y4    -----  #
  jdata.sim[, m1 := mean(y3), list(k,j2)]
  jdata.sim[, v1 := var(m1), list(k)]
  res$y3_Wa_Bk2 = jdata.sim[, mean(v1/v_y3_unc)]

  jdata.sim[, m1 := mean(y4), list(k,j2)]
  jdata.sim[, v1 := var(m1), list(k)]
  res$y4_Wa_Bk2 = jdata.sim[, mean(v1/v_y4_unc)]

  # -------  Wa,k' Bk of Y3,Y4 ------- #
  jdata.sim[, m1 := mean(y3), list(k,j1,j2)]
  jdata.sim[, v1 := var(m1), list(k,j2)]
  res$y3_Wak2_Bk1 = jdata.sim[, mean(v1/v_y3_unc)]

  jdata.sim[, m1 := mean(y4), list(k,j1,j2)]
  jdata.sim[, v1 := var(m1), list(k,j2)]
  res$y4_Wak2_Bk1 = jdata.sim[, mean(v1/v_y4_unc)]

  # ---------- Wa,k,k' of Y3,Y4 -------- #
  jdata.sim[, v1 := var(y3), list(k,j1,j2)]
  jdata.sim[is.na(v1),v1:=0]
  res$y3_Wak1k2 = jdata.sim[, mean(v1/v_y3_unc)]

  jdata.sim[, v1 := var(y4), list(k,j1,j2)]
  jdata.sim[is.na(v1),v1:=0]
  res$y4_Wak1k2 = jdata.sim[, mean(v1/v_y4_unc)]

  # ------- Wa Bk of Y3, Y4 -------- #
  jdata.sim[, m1 := mean(y3), list(k,j1)]
  jdata.sim[, v1 := var(m1), list(k)]
  res$y3_Wa_Bk1 = jdata.sim[, mean(v1/v_y3_unc)]

  jdata.sim[, m1 := mean(y4), list(k,j1)]
  jdata.sim[, v1 := var(m1), list(k)]
  res$y4_Wa_Bk1 = jdata.sim[, mean(v1/v_y4_unc)]

  # -------- network --------
  # first term
  jdata.sim[, m1 := mean(y3), list(k,j2)]
  jdata.sim[, a2 := mean(m1), list(k,j1)]
  jdata.sim[, v1 := var(a2), list(k)]
  res$y3net_Wa_Bk1 = jdata.sim[, mean(v1/v_y3_unc)]

  jdata.sim[, m1 := mean(y4), list(k,j2)]
  jdata.sim[, a2 := mean(m1), list(k,j1)]
  jdata.sim[, v1 := var(a2), list(k)]
  res$y4net_Wa_Bk1 = jdata.sim[, mean(v1/v_y4_unc)]

  res$y3var_diff = res$y3_Wa_Bk1- res$y3net_Wa_Bk1
  res$y4var_diff = res$y4_Wa_Bk1- res$y4net_Wa_Bk1

  return(data.frame(res))
}

generate_simualted_data = function() {
  load("inst/m2-mixt-y2003-main-fixb.rkiv")
  res_main = value
  model = res_main$model
  sim = m2.mixt.simulate.sim(model,fsize = 15)
  sdata = sim$sdata
  jdata = sim$jdata
  sdata$move=0
  jdata$move=1
  sdata[, birthyear := sample(1960:1980,.N,replace=T)]
  jdata[, birthyear := sample(1960:1980,.N,replace=T)]
  sdata[,wid:=sprintf("W%i",1:.N)]
  ns =sdata[,.N]
  jdata[,wid:=sprintf("W%i", ns+ (1:.N))]
  sdata[, ind1:= sample(c("Manufacturing","Services","Retail trade","Construction etc."),.N,replace=T)]
  jdata[, ind1:= sample(c("Manufacturing","Services","Retail trade","Construction etc."),.N,replace=T)]
  sdata[, educ:= sample(1:3,.N,replace=T)]
  jdata[, educ:= sample(1:3,.N,replace=T)]
  sdata[,size1 := .N,f1]
  jdata[,size1 := .N,f1]
  sdata[,va1 := exp(rnorm(.N)),f1]
  jdata[,va1 :=exp(rnorm(.N)),f1]
  save(sdata,jdata,file=sprintf("%s/data-tmp/tmp-2003-static.dat",local_opts$wdir))

  flog.info("!!! Using simulated data")
}
