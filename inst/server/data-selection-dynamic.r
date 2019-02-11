# this file creates the data set that will be used
# in the analysis of the paper
#
# related files are dataproc.r in the R folder in
# in the blm package.
require(assertthat)
require(blm)
require(foreign)
require(data.table)
require(SDMTools)

# load stata files, merge them
# loadDataRaw is a funciton in dataproc.R that extracts data 2 years before and after the anchoryear = 2003.
# In turn, loadDataRaw uses the proc loadData:
rdata = loadDataRaw(2003,path=sprintf("%s/data",local_opts$wdir),cache=F)
rdata_qtr=rdata[female==0]

# ====== AGGREGATE WORKER DATA TO YEARLY INFO ========
setkey(rdata_qtr,aret,wid) # order by year and worker ID
rdata_year  = rdata_qtr[, list(fids = length(unique(fid)),   # number of firms
                               m = sum(monthsworked),       # number of months worked
                               msd = sd(logrealmanadslon), # sd of earnings
                               fid = fid[1],
                               birthyear=birthyear[1], educ=educ[1],
                               y=logrealmanadslon[1], ind=industry[1],
                               va=valueadded[1], size = ant_anst[1] ,
                               no=.N),list(aret,wid)] # do it for each year and worker id,
assert_that(rdata_year[(m==12)&(fids)==1,max(msd)==0]) # fully employed worker should not have variation in earnings

# ======= SELECTING CONTINUING FIRMS ===========
# DEF: a continuing firm is a firm with a given worker fully employed in 2002 an 2004
# step1: select fully employed workers in both 2002 and 2004, in the same firm
setkey(rdata_year,wid,aret)
rdata_year_wide = rdata_year[(aret %in% c(2001,2002,2004,2005)) & (m==12) & (fids==1),list(f1=fid[1],f2=fid[2],f3=fid[3],f4=fid[4],.N,m=sum(m),aret1=aret[1],aret2=aret[2],aret3=aret[3],aret4=aret[4]),wid] # select the 2 years, extract month worked
assert_that(rdata_year_wide[,all(m<=48)])
assert_that(rdata_year_wide[N==4,all(aret1==2001)])
assert_that(rdata_year_wide[N==4,all(aret2==2002)])
assert_that(rdata_year_wide[N==4,all(aret3==2004)])
assert_that(rdata_year_wide[N==4,all(aret4==2005)])
fids_continuing = rdata_year_wide[(f1==f2)&(f1==f3)&(f1==f4)&(m==48),unique(f1)] # select firms with at least one stayers, fully employed in both periods

# drop other firms
rdata_qtr  = rdata_qtr[ fid %in% fids_continuing]
rdata_year = rdata_year[fid %in% fids_continuing]
assert_that(all(rdata_qtr[,sort(unique(fid))]==rdata_year[,sort(unique(fid))]))
assert_that(rdata_qtr[,length(unique(wid))]>=rdata_year[,length(unique(wid))])

# ======= SELECTING FULLY EMPLOYED WORKERS ===========
# keep only m==12 and fids==1 in 2001/2 and on 2004/5,
setkey(rdata_year,wid,aret)
rdata_year_wide = rdata_year[(aret %in% c(2001,2002,2004,2005)) & (m==12) & (fids==1),list(f1=fid[1],f2=fid[2],f3=fid[3],f4=fid[4],.N,m=sum(m),aret1=aret[1],aret2=aret[2],aret3=aret[3],aret4=aret[4]),wid] # select the 2 years, extract month worked
wids_fullemployed = rdata_year_wide[(f1==f2)&(f3==f4)&(m==48),wid]

rdata_qtr  = rdata_qtr [aret %in% c(2001,2002,2003,2004,2005)][wid %in% wids_fullemployed]
rdata_year = rdata_year[aret %in% c(2001,2002,2003,2004,2005)][wid %in% wids_fullemployed]

# ====== CONSTRUCT INFO ABOUT WHAT HAPPENS DURING 2003 ========

setkey(rdata_qtr,wid,aret,quarter)
rdata_qtr[,fid.l1q4 := rdata_qtr[J(wid,aret-1,4),fid]] # get fid in previous year quarter 1
rdata_qtr[,fid.f1q1 := rdata_qtr[J(wid,aret+1,1),fid]] # get fid in following year quarter 4
rdata_qtr[,fid.q1   := rdata_qtr[J(wid,aret,1),fid]]   # get fid in current year quarter 1
rdata_qtr[,fid.q2   := rdata_qtr[J(wid,aret,2),fid]]   # get fid in current year quarter 2
rdata_qtr[,fid.q3   := rdata_qtr[J(wid,aret,3),fid]]   # get fid in current year quarter 3
rdata_qtr[,fid.q4   := rdata_qtr[J(wid,aret,4),fid]]   # get fid in current year quarter 4

assert_that(rdata_qtr[aret==2003,all(!is.na(fid.l1q4))])
assert_that(rdata_qtr[aret==2003,all(!is.na(fid.f1q1))])

# Individual is the unit of observation
mdata2003 = rdata_qtr[(aret==2003),list(u=4-sum(employed),
                                        fcount =length(unique(fid)), # number of firm identifiers in 2003
                                        fcount2=length(unique(c(fid,fid.l1q4,fid.f1q1))),   #number of firm identifiers between 2002,2003,2004 conditional on being full year employed
                                        m00=(fid.q1[1]==fid.f1q1[1]) & (fid.q2[1]==fid.f1q1[1]) & (fid.q3[1]==fid.f1q1[1]) & (fid.q4[1]==fid.f1q1[1]) & (fid.q1[1]==fid.l1q4[1])  , # (f1,)f2,f2,f2,f2(,f2)
                                        m01=(fid.q1[1]==fid.f1q1[1]) & (fid.q2[1]==fid.f1q1[1]) & (fid.q3[1]==fid.f1q1[1]) & (fid.q4[1]==fid.f1q1[1]) & (fid.q1[1]!=fid.l1q4[1])  , # (f1,)f2,f2,f2,f2(,f2)
                                        m12=(fid.q1[1]==fid.l1q4[1]) & (fid.q2[1]==fid.f1q1[1]) & (fid.q3[1]==fid.f1q1[1]) & (fid.q4[1]==fid.f1q1[1]) & (fid.f1q1[1]!=fid.l1q4[1]), # (f1,)f1,f2,f2,f2(,f2)
                                        m23=(fid.q1[1]==fid.l1q4[1]) & (fid.q2[1]==fid.l1q4[1]) & (fid.q3[1]==fid.f1q1[1]) & (fid.q4[1]==fid.f1q1[1]) & (fid.f1q1[1]!=fid.l1q4[1]), # (f1,)f1,f1,f2,f2(,f2)
                                        m34=(fid.q1[1]==fid.l1q4[1]) & (fid.q2[1]==fid.l1q4[1]) & (fid.q3[1]==fid.l1q4[1]) & (fid.q4[1]==fid.f1q1[1]) & (fid.f1q1[1]!=fid.l1q4[1]), # (f1,)f1,f1,f1,f2(,f2)
                                        m45=(fid.q1[1]==fid.l1q4[1]) & (fid.q2[1]==fid.l1q4[1]) & (fid.q3[1]==fid.l1q4[1]) & (fid.q4[1]==fid.l1q4[1]) & (fid.q4[1]!=fid.f1q1[1])  , # (f1,)f1,f1,f1,f1(,f2)                                    move=(fid.l1q4[1]!=fid.f1q1[1]),
                                        move=(fid.f1q1[1]!=fid.l1q4[1]),
                                        present_in_2003 = TRUE,
                                        qtr_in_2003=.N),wid]
mdata2003[qtr_in_2003<4,c('m00','m01','m12','m23','m34','m45'):=FALSE,with=FALSE]

# add workers that have 0 quarters in 2003 and full years in 2002/2004
# Identified as the difference between workers fully employed in 2002 and 2004 that are present in 2003.
wids_noinfo = setdiff(rdata_qtr[,unique(wid)],rdata_qtr[aret %in% c(2003),unique(wid)])

# Individual worker is the unit of observation:
mdata2003_noinfo = rdata_qtr[wid %in% wids_noinfo,list(u=0,
                                                       fcount =0,
                                                       fcount2=length(unique(fid)),
                                                       m01=FALSE,m12=FALSE,m23=FALSE,m34=FALSE,m45=FALSE,m00=FALSE,
                                                       move=(fid[(aret==2002)&(quarter==4)]!=fid[(aret==2004)&(quarter==1)]),
                                                       present_in_2003 = FALSE,
                                                       qtr_in_2003=0),wid]
mdata2003 = rBind(mdata2003,mdata2003_noinfo)

assert_that(mdata2003[,all(is.na(move)==FALSE)])
assert_that(mdata2003[m01|m12|m23|m34|m45,all(qtr_in_2003==4)])
assert_that(mdata2003[m01|m12|m23|m34|m45,all(m01+m12+m23+m34+m45<=1)])
assert_that(mdata2003[m01|m12|m23|m34|m45,all(move==TRUE)])

# remove 2003, we have collected all necessary info
rdata_year = rdata_year[aret %in% c(2001,2002,2004,2005)]

# ----- construct movers -------
setkey(rdata_year, wid, aret)
widm = rdata_year[,list(fid[2]!=fid[3]),wid][V1==TRUE,wid]
assert_that(length(widm) == mdata2003[,sum(move==TRUE)])
assert_that(all( sort(widm) == mdata2003[move==TRUE,sort(unique(wid))])) # This checks that the two ways of constructing the movers coincide.

# make a wide table
rdatam = rdata_year[wid %in% widm]
setkey(rdatam,wid,aret)
jdata = rdatam[, list(y1=y[1],aret1=aret[1],
                      y2=y[2],f1=fid[2],aret2=aret[2],
                      y3=y[3],f2=fid[3],aret3=aret[3],
                      y4=y[4],aret4=aret[4],
                      birthyear=birthyear[1],educ=educ[1],
                      ind1=ind[2],size1=size[2],va1=va[2],
                      ind2=ind[3],size2=size[3],va2=va[3]),wid]

assert_that(jdata[,sum(is.na(f2))==0])
assert_that(jdata[,sum(is.na(f1))==0])
assert_that(jdata[,sum(is.na(y1))==0])
assert_that(jdata[,sum(is.na(y2))==0])
assert_that(jdata[,sum(is.na(y3))==0])
assert_that(jdata[,sum(is.na(y4))==0])
assert_that(jdata[,all(aret1==2001)])
assert_that(jdata[,all(aret2==2002)])
assert_that(jdata[,all(aret3==2004)])
assert_that(jdata[,all(aret4==2005)])

# merge in 2003 info
setkey(mdata2003,wid)
setkey(jdata,wid)
jdata = mdata2003[jdata]

# ----- stayers & movers ---

# make a wide table
setkey(rdata_year,wid,aret)
sdata = rdata_year[, list(y1=y[1],aret1=aret[1],
                          y2=y[2],f1=fid[2],aret2=aret[2],
                          y3=y[3],f2=fid[3],aret3=aret[3],
                          y4=y[4],aret4=aret[4],
                          birthyear=birthyear[1],educ=educ[1],
                          ind1=ind[2],size1=size[2],va1=va[2],
                          ind2=ind[3],size2=size[3],va2=va[3]),wid]

# merge in the moving info
setkey(mdata2003,wid)
setkey(sdata,wid)
sdata = mdata2003[sdata]


# for each data set we compute the following statistics
get.stats <- function(data,movers=FALSE) {
  rr = list()
  rr$nwid       = data[,length(unique(wid))]   # unique worker
  rr$nfid       = data[,length(unique(f1))]

  # creating firm info from 2002
  fdata = data[,list(.N,ind=ind1[1],size1=size1[1],va1=va1[1],Nm=sum(move==TRUE)),f1]

  data[,asize := .N,f1]

  rr$nfirm_actualsize_ge10       = fdata[N>=10,.N]
  rr$nfirm_actualsize_ge50       = fdata[N>=50,.N]
  rr$nfirm_reportedsize_ge10       = fdata[size1>=10,.N]
  rr$nfirm_reportedsize_ge50       = fdata[size1>=50,.N]
  rr$nfirm_movers_ge1       = fdata[Nm>=1,.N]
  rr$nfirm_movers_ge5       = fdata[Nm>=5,.N]
  rr$nfirm_movers_ge10      = fdata[Nm>=10,.N]
  rr$firm_reportedsize_mean   = fdata[,mean(size1)]
  rr$firm_reportedsize_median = fdata[,median(size1)]
  rr$firm_actualsize_mean   = fdata[,mean(N)]
  rr$firm_actualsize_median = fdata[,median(N)]
  rr$firm_reportedsize_median_worker = data[,median(size1)]
  rr$firm_actualsize_median_worker = data[,median(asize)]
  rr$worker_share_ind_manu   = data[,mean(ind1=="Manufacturing")]
  rr$worker_share_ind_serv   = data[,mean(ind1=="Services")]
  rr$worker_share_ind_retail   = data[,mean(ind1=="Retail trade")]
  rr$worker_share_ind_cons   = data[,mean(ind1=="Construction etc.")]
  rr$firm_mean_log_va = fdata[va1>0,mean(log(va1))]
  rr$firm_var_log_va = fdata[va1>0,var(log(va1))]
  rr$firm_neg_va = fdata[va1<=0,.N]

  rr$worker_share_educ1   = data[,mean(educ==1)]
  rr$worker_share_educ2   = data[,mean(educ==2)]
  rr$worker_share_educ3   = data[,mean(educ==3)]
  rr$worker_var_log_wage  = data[,var(y2)]
  rr$worker_mean_log_wage = data[,mean (y2)]

  # between firm variance
  data[,fmw := mean(y2),f1]
  rr$between_firm_wage_var = data[, var(fmw)]

  data[,age:=2002 - birthyear]
  rr$worker_share_age_0_30   = data[,mean(age<=30)]
  rr$worker_share_age_31_50  = data[,mean(age>=31 & age<=50)]
  rr$worker_share_age_51_inf = data[,mean(age>=51)]
  rr$worker_stay_all = data[,sum((move==FALSE)&(m00==TRUE))]

  return(rr)
}

rrs = get.stats(sdata)
rrj = get.stats(jdata,movers=TRUE)

# compute between class variance of earnings and va
sdata[,list(mean(y1),.N),j1][,wt.var(V1,N)]


rr = jdata[m01|m12|m23|m34|m45,list(y1=mean(y1),y2=mean(y2),y3=mean(y3),y4=mean(y4),.N,m=m01+2*m12+3*m23+4*m34+5*m45),list(m01,m12,m23,m34,m45,fcount2)][order(m)]
rr = melt(rr,id.vars = c("m",'N'),measure.vars = c('y1','y2','y3','y4'))
ggplot(rr,aes(x=m,y=value,color=factor(variable))) + geom_line() +geom_point(aes(size=N))+theme_bw()
ggplot(rr,aes(x=factor(variable),y=value,color=factor(m),group=factor(m))) + geom_line() +geom_point(aes(size=N))+theme_bw()

save(rrs,rrj,file=sprintf("%s/data-tmp/tab-summary-dynamic.dat",local_opts$wdir))
save(sdata,jdata,file=sprintf("%s/data-tmp/data-dynamic.dat",local_opts$wdir))





