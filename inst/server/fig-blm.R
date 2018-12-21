# This file generates all the figures of the paper and the
# appendix.

fig.main <- function() {

  # ======= PAPER ========= #
  fig.static.mixt.means() # FIGURE 2

  # ======= SUPPLMENT ======== #
  plot.static.proba.gibbs() # Figure S7
  fig.dynamic.mixt.means()  # Figure S14
  fig.dynamic.mixture.fit()

  # ======= RESPONSE======== #
  figure.hybrid()   # Figure E9


  # shimer-smith
  fig.shimersmith.wages()

  # S16
  # SS17
}

tables.all <- function() {

  # ======= PAPER ========= #
  table.static.clusters()       # Table 1
  table.static.mixt.vdec()      # Table 2
  table.m4.parameters()         # Table 3
  table.endogeneousMobility()   # Table 4
  table.statedependence()       # Table 5

  # ======= SUPPLEMENT ========= #
  table.static.mixt.vdec()        # Table S10
  table.dynamic.mixt.vdec()       # Table S16, S17
  table.m2.mixt.vdec.robust.nk()  # Table S4
  table.dynamic.clusters()                # Table S11


  # Reclassification table
  tab.dirty_iteration()



  # Ednogenous mobility

  # ==== RESPONSE ==== #
  table.cluster.extra.moments()  # Table 1
  fig.shimersmith.narrow()       # table 5 Shimer Smith

}






# ======================== DATA STATS =====================================

#' this generates the summary statistics by cluster
#' which ends up being table 2 in the paper
table.static.clusters <- function() {
  # create table for static
  m2stats = res.load("m2-stats-d2003")
  gstats = m2stats$gstats

  tt =
    tt_rule_top() +
    tt_text_row(c("class:",paste(1:10),"all")) +
    tt_rule_mid() +
    tt_text_row("number of workers") %:% tt_numeric_row(gstats$nwid,dec=0) +
    tt_text_row("number of firms") %:% tt_numeric_row(gstats$nfid,dec=0) +
    tt_text_row("mean firm reported size               ")  %:% tt_numeric_row(gstats$firm_reportedsize_mean,dec=2) +
    tt_text_row("number of firms $\\geq10$ (actual size)") %:% tt_numeric_row(gstats$nfirm_actualsize_ge10,dec=0) +
    tt_text_row("number of firms $\\geq50$ (actual size)") %:% tt_numeric_row(gstats$nfirm_actualsize_ge50,dec=0) + tt_spacer_row(9) +
    tt_text_row("\\% high school drop out               ") %:% tt_numeric_row(100*gstats$worker_share_educ1,perc=T,dec=1) +
    tt_text_row("\\% high school graduates              ") %:% tt_numeric_row(100*gstats$worker_share_educ2,perc=T,dec=1) +
    tt_text_row("\\% some college                       ") %:% tt_numeric_row(100*gstats$worker_share_educ3,perc=T,dec=1) + tt_spacer_row(9) +
    tt_text_row("\\% workers younger than 30            ") %:% tt_numeric_row(100*gstats$worker_share_age_0_30,perc=T,dec=1) +
    tt_text_row("\\% workers between 31 and 50          ") %:% tt_numeric_row(100*gstats$worker_share_age_31_50,perc=T,dec=1) +
    tt_text_row("\\% workers older than 51              ") %:% tt_numeric_row(100*gstats$worker_share_age_51_inf,perc=T,dec=1) + tt_spacer_row(9) +
    tt_text_row("\\% workers in manufacturing           ") %:% tt_numeric_row(100*gstats$worker_share_ind_manu,perc=T,dec=1) +
    tt_text_row("\\% workers in services                ") %:% tt_numeric_row(100*gstats$worker_share_ind_serv,perc=T,dec=1) +
    tt_text_row("\\% workers in retail and trade        ") %:% tt_numeric_row(100*gstats$worker_share_ind_retail,perc=T,dec=1) +
    tt_text_row("\\% workers in construction            ") %:% tt_numeric_row(100*gstats$worker_share_ind_cons,perc=T,dec=1) + tt_spacer_row(9) +
    tt_text_row("mean log-earnings                     ")  %:% tt_numeric_row(gstats$worker_mean_log_wage,dec=2) +
    tt_text_row("variance of log-earnings              ")  %:% tt_numeric_row(gstats$worker_var_log_wage,dec=3) +
    tt_text_row("skewness of log-earnings              ")  %:% tt_numeric_row(gstats$worker_skewness_log_wage,dec=3) +
    tt_text_row("kurtosis of log-earnings              ")  %:% tt_numeric_row(gstats$worker_kurtosis_log_wage,dec=3) +
    tt_text_row("between-firm variance of log-earnings ")  %:% tt_numeric_row(gstats$between_firm_wage_var,dec=4) +
    tt_text_row("mean log-value-added per worker       ")  %:% tt_numeric_row(gstats$firm_mean_log_va,dec=2) +
    tt_rule_bottom()

  tab = tt_tabularize(tt,"l rrrrrrrrrr|r")
  tt_save(tab,filename="inst/figures/tab-summary-clusters-m2-tmp.tex",stand_alone=F)

  brew("inst/templates/tab-summary-clusters.br","inst/fig/build/tab-summary-clusters-m2.tex")
  file.copy("inst/fig/build/tab-summary-clusters-m2.tex","inst/figures",overwrite = T)




  # create table for dynamic
  m4stats = res.load("m4-stats-d2003")
  gstats  = m4stats$gstats
  brew("inst/templates/tab-summary-clusters.br","inst/fig/build/tab-summary-clusters-m4.tex")
  file.copy("inst/fig/build/tab-summary-clusters-m4.tex","inst/figures",overwrite=T)
}

table.dynamic.clusters <- function() {
  # create table for static
  m2stats = res.load("m4-stats-d2003")
  gstats = m2stats$gstats

  tt =
    tt_rule_top() +
    tt_text_row(c("class:",paste(1:10),"all")) +
    tt_rule_mid() +
    tt_text_row("number of workers") %:% tt_numeric_row(gstats$nwid,dec=0) +
    tt_text_row("number of firms") %:% tt_numeric_row(gstats$nfid,dec=0) +
    tt_text_row("mean firm reported size               ")  %:% tt_numeric_row(gstats$firm_reportedsize_mean,dec=2) +
    tt_text_row("number of firms $\\geq10$ (actual size)") %:% tt_numeric_row(gstats$nfirm_actualsize_ge10,dec=0) +
    tt_text_row("number of firms $\\geq50$ (actual size)") %:% tt_numeric_row(gstats$nfirm_actualsize_ge50,dec=0) + tt_spacer_row(9) +
    tt_text_row("\\% high school drop out               ") %:% tt_numeric_row(100*gstats$worker_share_educ1,perc=T,dec=1) +
    tt_text_row("\\% high school graduates              ") %:% tt_numeric_row(100*gstats$worker_share_educ2,perc=T,dec=1) +
    tt_text_row("\\% some college                       ") %:% tt_numeric_row(100*gstats$worker_share_educ3,perc=T,dec=1) + tt_spacer_row(9) +
    tt_text_row("\\% workers younger than 30            ") %:% tt_numeric_row(100*gstats$worker_share_age_0_30,perc=T,dec=1) +
    tt_text_row("\\% workers between 31 and 50          ") %:% tt_numeric_row(100*gstats$worker_share_age_31_50,perc=T,dec=1) +
    tt_text_row("\\% workers older than 51              ") %:% tt_numeric_row(100*gstats$worker_share_age_51_inf,perc=T,dec=1) + tt_spacer_row(9) +
    tt_text_row("\\% workers in manufacturing           ") %:% tt_numeric_row(100*gstats$worker_share_ind_manu,perc=T,dec=1) +
    tt_text_row("\\% workers in services                ") %:% tt_numeric_row(100*gstats$worker_share_ind_serv,perc=T,dec=1) +
    tt_text_row("\\% workers in retail and trade        ") %:% tt_numeric_row(100*gstats$worker_share_ind_retail,perc=T,dec=1) +
    tt_text_row("\\% workers in construction            ") %:% tt_numeric_row(100*gstats$worker_share_ind_cons,perc=T,dec=1) + tt_spacer_row(9) +
    tt_text_row("mean log-earnings                     ")  %:% tt_numeric_row(gstats$worker_mean_log_wage,dec=2) +
    tt_text_row("variance of log-earnings              ")  %:% tt_numeric_row(gstats$worker_var_log_wage,dec=3) +
    tt_text_row("skewness of log-earnings              ")  %:% tt_numeric_row(gstats$worker_skewness_log_wage,dec=3) +
    tt_text_row("kurtosis of log-earnings              ")  %:% tt_numeric_row(gstats$worker_kurtosis_log_wage,dec=3) +
    tt_text_row("between-firm variance of log-earnings ")  %:% tt_numeric_row(gstats$between_firm_wage_var,dec=4) +
    tt_text_row("mean log-value-added per worker       ")  %:% tt_numeric_row(gstats$firm_mean_log_va,dec=2) +
    tt_rule_bottom()

  tab = tt_tabularize(tt,"l rrrrrrrrrr|r")
  tt_save(tab,filename="inst/figures/tab-summary-clusters-m4.tex",stand_alone=F)

  brew("inst/templates/tab-summary-clusters.br","inst/fig/build/tab-summary-clusters-m2.tex")
  file.copy("inst/fig/build/tab-summary-clusters-m2.tex","inst/figures",overwrite = T)




  # create table for dynamic
  m4stats = res.load("m4-stats-d2003")
  gstats  = m4stats$gstats
  brew("inst/templates/tab-summary-clusters.br","inst/fig/build/tab-summary-clusters-m4.tex")
  file.copy("inst/fig/build/tab-summary-clusters-m4.tex","inst/figures",overwrite=T)
}




table.cluster.extra.moments <- function() {
  # create table with extra moments
  m2stats = res.load("m2-stats-d2003")
  gstats = m2stats$gstats[1:10,]

  tt = tt_text_row(c("cluster","N","mean","var.","skew.","kurt.")) +
  tt_rule_mid_partial(list(c(2,6))) +
  tt_numeric_column(gstats$j1,dec=0) %&%
      tt_numeric_column(gstats$nwid,dec=0)  %&%
      tt_numeric_column(gstats$worker_mean_log_wage,dec=2)  %&%
      tt_numeric_column(gstats$worker_var_log_wage,dec=2) %&%
      tt_numeric_column(gstats$worker_skewness_log_wage,dec=2) %&%
      tt_numeric_column(gstats$worker_kurtosis_log_wage,dec=1)

  tab = tt_tabularize(tt,header = "rrrrrr", pretty_rules=T)
  tt_save(tab,filename="inst/figures/tab-cluster-extra-moments.tex",stand_alone=F)

}


#' creates the table with mover counts
table.movers.count <- function() {
  m2stats = res.load("m2-stats-d2003")
  cstats = m2stats$cstats
  mstats = m2stats$mstats
  brew('inst/templates/tab-summary-mobility.br','inst/fig/build/tab-summary-mobility.tex')
  file.copy('inst/fig/build/tab-summary-mobility.tex',"inst/figures",overwrite = T)

  m2stats = res.load("m4-stats-d2003")
  cstats = m2stats$cstats
  mstats = m2stats$mstats
  brew('inst/templates/tab-summary-mobility.br','inst/fig/build/tab-summary-mobility-m4.tex')
  file.copy('inst/fig/build/tab-summary-mobility-m4.tex',"inst/figures",overwrite = T)
}

table.movers.wages <- function() {
  mstats = res.load("m2-stats-d2003")$mstats
  #colnames(mstats) = c("j1","j2","m1","sd1","m2","sd2","m3","sd3","m4","sd4","N")

  mstats= data.table(mstats)
  mstats[,j1b := ceiling(j1/2)]
  rr1 = mstats[,wt.mean(m1,N),list(j1b,j2)]

  ggplot(rr1,aes(x=factor(j2),y=V1,group=j1b,color=factor(j1b))) +
    geom_line() + theme_bw() + xlab("firm class k in period 2") + ylab("mean log earnings") +
    theme(legend.position = "none") +  coord_cartesian(ylim=c(9.6,10.75))
  ggsave("inst/figuresfig-static-k2_on_y1.pdf",dpi=100,width = 6,height = 5)


  mstats[,j2b := ceiling(j2/2)]
  rr2 = mstats[,wt.mean(m2,N),list(j1,j2b)]

  ggplot(rr2,aes(x=factor(j1),y=V1,group=j2b,color=factor(j2b))) +
    geom_line() + theme_bw() + xlab("firm class k in period 1") + ylab("mean log earnings") +
    theme(legend.position = "none") +  coord_cartesian(ylim=c(9.6,10.75))
  ggsave("inst/figuresfig-static-k1_on_y2.pdf",dpi=100,width = 6,height = 5)
}

# =========================  STATIC -- MIXTURE =====================

#' plotting the means and proportions for the main estimate
#' of the mixture model.
fig.static.mixt.means <- function() {
  res_main    = res.load("m2-mixt-y2003-main-fixb")
  res_bs      = res.load("m2-mixt-d2003-bootstrap")$mixt_all
  cstats      = res.load("m2-stats-d2003")$cstats
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  res_mixt_bs = data.table(ldply(model_mixt_bs,function(x) {
    I  = order(colSums(-x$model$A1))
    If = order(rowSums(x$model$A1*x$model$pk0[1,,]))
    melt(x$model$A1[,I][If,],c('l','k')) }
  ))
  res_sd      = res_mixt_bs[,list(q0=quantile(value,0.025),q1=quantile(value,0.975),m=mean(value)),list(l,k)]

  I = order(colSums(-res_main$model$A1))
  res_sd[,value:=res_main$model$A1[l,I[k]] ,list(k,l)]

  lmin = min(cstats[,m1-3*sd1])
  lmax = max(cstats[,m1+3*sd1])
  cstats$k=1

  gp = wplot(res_main$model$A1[,I]) + xlab("firm class k") +theme_bw() +
    scale_y_continuous("log-earnings")+theme(legend.position = "none") +
    geom_errorbar(aes(ymin=q0,ymax=q1),data=res_sd,width=0.2) +
    scale_color_brewer(palette="Spectral")

  ggsave("../figures/fig-mixt2-d2003-wage-uc.pdf",gp,width = 6.5,height = 5)

  # proportions
  dpk1 = m2.get.pk1(res_main$model)
  pk_m  = acast(dpk1[,sum(pr_j1j2k),list(j1,k)],j1~k,value.var = "V1")
  NNs = res_main$model$NNs*10 # used 10% sample
  NNm = res_main$model$NNm
  share_s  = sum(NNs)/(sum(NNm) + sum(NNs))
  pk_unc = share_s*rdim(res_main$model$pk0[,,I],res_main$model$nf,res_main$model$nk) +
             (1- share_s) * pk_m

  dd = melt(pk_unc,c('l','k'))
  gp = ggplot(dd,aes(x=factor(l),y=value,fill=factor(k))) +
    geom_bar(position="stack",stat = "identity") + theme_bw() +
    theme(legend.position = "none") +
    xlab("firm class k") + ylab("type proportions") +
    scale_fill_brewer(palette="Spectral")

  ggsave("../figures/fig-mixt2-d2003-pk.pdf",gp,width = 6.5,height = 5)
}

table.static.mixt.vdec <- function() {

  res_main       = res.load("m2-mixt-y2003-main-fixb")
  model_mixt_bs  = m2.getboostrap()$mixt_all
  # add additional boostraps

  vdec.extract <- function(vdec) {
    VV = vdec$cc
    r = data.frame(var_y=VV[1,1],var_k=VV[2,2],var_l=VV[3,3],var_e=VV[4,4],cov_kl=VV[2,3])
    r$cor_kl        = r$cov_kl/sqrt(r$var_k*r$var_l)
    r$var_l_tshare  = r$var_l/r$var_y
    r$var_k_tshare  = r$var_k/r$var_y
    r$var_e_tshare  = r$var_e/r$var_y
    r$cov_kl_tshare = 2*r$cov_kl/r$var_y
    r
  }

  df2list <- function(dd,values,names="variable") {
    res = as.list(dd[,get(values)])
    names(res) = dd[,get(names)]
    res
  }

  # compute boostrap s.e.
  res_mixt_bs = data.table(ldply(model_mixt_bs,function(x) vdec.extract(x$vdec)))
  res_mixt_bs = data.table(melt(res_mixt_bs,id=".id"))
  res_sd      = res_mixt_bs[,list(q0=quantile(value,0.025),q1=quantile(value,0.975),m=mean(value),sd=sd(value)),variable]
  main_sd     = df2list(res_sd,"sd")
  bs_vdec_m1  = df2list(res_sd,"m")
  bs_vdec_q0  = df2list(res_sd,"q0")
  bs_vdec_q1  = df2list(res_sd,"q1")

  main = vdec.extract(res_main$vdec)

  # ===== EXTRACT REALLOCATION RESULTS ====== #
  rr_mean_effects = res.load("m2-mixt-d2003-meffects-leg2")
  res_means         = df2list(rr_mean_effects$diff,"main")
  res_means_sd      = df2list(rr_mean_effects$diff,"sd")
  res_means_m1      = df2list(rr_mean_effects$diff,"m1")
  res_means_q0      = df2list(rr_mean_effects$diff,"q0")
  res_means_q1      = df2list(rr_mean_effects$diff,"q1")

  require(textables)

  tt =
    tt_rule_top() +
    tt_text_row("{ \\bf Variance decomposition ($\\times 100$) }",cspan = 5,center = "c") + tt_spacer_row(-3) +
	  tt_text_row(c("$\\frac{Var(\\alpha)}{Var(y)}$",
	  "$\\frac{Var(\\psi)}{Var(y)}$",
	  "$\\frac{2 Cov(\\alpha,\\psi)}{Var(y)}$",
	  "$\\frac{Var(\\varepsilon)}{Var(y)}$",
	  "$Corr(\\alpha,\\psi)$")) +
    tt_rule_mid() +
  tt_numeric_row(
    c(100*main$var_k_tshare ,
      100*main$var_l_tshare,
      100*main$cov_kl_tshare,
      100*main$var_e_tshare,
      100*main$cor_kl),
    percent=F,dec=2) +
  tt_spacer_row(-9) +
  tt_numeric_row(
	    c(100*main_sd$var_k_tshare,
	      100*main_sd$var_l_tshare,
	      100*main_sd$cov_kl_tshare,
	      100*main_sd$var_e_tshare,
	      100*main_sd$cor_kl),
	    percent=F,dec=2,surround="{\\scriptsize (%s)}")  + tt_spacer_row(6) +
	tt_text_row(" { \\bf Reallocation exercise ($\\times 100$) }",cspan = 5,center = "c") + tt_spacer_row(-5) +
	  tt_text_row( c("Mean","Median","10\\%-quantile", "90\\%-quantile", "Variance")) + tt_spacer_row(-5) +
    tt_rule_mid() +
	  tt_numeric_row(100*as.numeric(res_means[c("mean","median","d1","d9","var")]),dec=2) +
	  tt_spacer_row(-9) +
    tt_numeric_row(100*as.numeric(res_means_sd[c("mean","median","d1","d9","var")]),dec=2,surround="{\\scriptsize (%s)}") +
    tt_rule_bottom()

	tab = tt_tabularize(tt,header = "ccccc", pretty_rules=F)
	tt_save(tab,filename="inst/figures/tab-static-main.tex",stand_alone=F)
	tt_save(tab,filename="~/tmp/tmp.tex",stand_alone=T)

	vnames = c("var_k_tshare","var_l_tshare","cov_kl_tshare","var_e_tshare","cor_kl")

	# SECOND TABLE
	tt =
	  tt_rule_top() +
	  TR("") %:% tt_text_row("{ \\bf Variance decomposition ($\\times 100$) }",cspan = 5,center = "c") + tt_spacer_row(-3) +
	  tt_text_row(c("","$\\frac{Var(\\alpha)}{Var(y)}$",
	                "$\\frac{Var(\\psi)}{Var(y)}$",
	                "$\\frac{2 Cov(\\alpha,\\psi)}{Var(y)}$",
	                "$\\frac{Var(\\varepsilon)}{Var(y)}$",
	                "$Corr(\\alpha,\\psi)$")) +
	  tt_rule_mid() +
	  TR("estimate") %:% tt_numeric_row(100*as.numeric(main[vnames]),percent=F,dec=2) +
	  tt_rule_mid_partial(list(c(2,6))) +
	  TR("bootstrap mean")             %:% TR(100*as.numeric(bs_vdec_m1[vnames]),percent=F,dec=1)  +
	  TR("bootstrap 2.5\\%-quantile")  %:% TR(100*as.numeric(bs_vdec_q0[vnames]),percent=F,dec=1)  +
	  TR("bootstrap 97.5\\%-quantile") %:% TR(100*as.numeric(bs_vdec_q1[vnames]),percent=F,dec=1)  +
	  tt_spacer_row(6) +
	  TR("") %:% tt_text_row(" { \\bf Reallocation exercise ($\\times 100$) }",cspan = 5,center = "c") + tt_spacer_row(-5) +
	  TR("") %:% tt_text_row( c("Mean","Median","10\\%-quantile", "90\\%-quantile", "Variance")) + tt_spacer_row(-5) +
	  tt_rule_mid() +
	  TR("estimate") %:% TR(100*as.numeric(res_means[c("mean","median","d1","d9","var")]),dec=2)  +
	  tt_rule_mid_partial(list(c(2,6))) +
	  TR("bootstrap mean")             %:% TR(100*as.numeric(res_means_m1[c("mean","median","d1","d9","var")]),dec=2)  +
	  TR("bootstrap 2.5\\%-quantile")  %:% TR(100*as.numeric(res_means_q0[c("mean","median","d1","d9","var")]),dec=2)  +
	  TR("bootstrap 97.5\\%-quantile") %:% TR(100*as.numeric(res_means_q1[c("mean","median","d1","d9","var")]),dec=2)  +
	tt_rule_bottom()

	tab = tt_tabularize(tt,header = "lccccc", pretty_rules=F)
	tt_save(tab,filename="inst/figures/tab-static-main-bootstrap-details.tex",stand_alone=F)

}

table.dynamic.mixt.vdec <- function() {

  res_main        = res.load("m4-mixt-d2003-main")
  model_mixt_bs   = m4.getboostrap()

  vdec.extract <- function(vdec) {
    VV = vdec$cc
    r = data.frame(var_y=VV[1,1],var_k=VV[2,2],var_l=VV[3,3],var_e=VV[4,4],cov_kl=VV[2,3])
    r$cor_kl        = r$cov_kl/sqrt(r$var_k*r$var_l)
    r$var_l_tshare  = r$var_l/r$var_y
    r$var_k_tshare  = r$var_k/r$var_y
    r$var_e_tshare  = r$var_e/r$var_y
    r$cov_kl_tshare = 2*r$cov_kl/r$var_y
    r
  }

  df2list <- function(dd,values,names="variable") {
    res = as.list(dd[,get(values)])
    names(res) = dd[,get(names)]
    res
  }

  # compute boostrap s.e.
  res_mixt_bs = data.table(ldply(model_mixt_bs$mixt_all,function(x) vdec.extract(x$vdec)))
  res_mixt_bs = data.table(melt(res_mixt_bs,id=".id"))
  res_sd      = res_mixt_bs[,list(q0=quantile(value,0.025),q1=quantile(value,0.975),m=mean(value),sd=sd(value)),variable]
  main_sd     = df2list(res_sd,"sd")
  bs_vdec_m1  = df2list(res_sd,"m")
  bs_vdec_q0  = df2list(res_sd,"q0")
  bs_vdec_q1  = df2list(res_sd,"q1")

  main = vdec.extract(res_main$vdec)

  # ===== EXTRACT REALLOCATION RESULTS ====== #
  rr_mean_effects = res.load("m4-mixt-d2003-meffects-leg2")
  res_means         = df2list(rr_mean_effects$diff,"main")
  res_means_sd      = df2list(rr_mean_effects$diff,"sd")
  res_means_m1      = df2list(rr_mean_effects$diff,"m1")
  res_means_q0      = df2list(rr_mean_effects$diff,"q0")
  res_means_q1      = df2list(rr_mean_effects$diff,"q1")

  require(textables)
  vnames = c("var_k_tshare","var_l_tshare","cov_kl_tshare","var_e_tshare","cor_kl")

  tt =
    tt_rule_top() +
    tt_text_row("{ \\bf Variance decomposition ($\\times 100$) }",cspan = 5,center = "c") + tt_spacer_row(-3) +
    tt_text_row(c("$\\frac{Var(\\alpha)}{Var(y)}$",
                  "$\\frac{Var(\\psi)}{Var(y)}$",
                  "$\\frac{2 Cov(\\alpha,\\psi)}{Var(y)}$",
                  "$\\frac{Var(\\varepsilon)}{Var(y)}$",
                  "$Corr(\\alpha,\\psi)$")) +
    tt_rule_mid() +
    tt_numeric_row(100*as.numeric(main[vnames]),percent=F,dec=2) +
    tt_spacer_row(-9) +
    tt_numeric_row(100*as.numeric(main_sd[vnames]),percent=F,dec=2,surround="{\\scriptsize (%s)}")  + tt_spacer_row(6) +
    tt_text_row(" { \\bf Reallocation exercise ($\\times 100$) }",cspan = 5,center = "c") + tt_spacer_row(-5) +
    tt_text_row( c("Mean","Median","10\\%-quantile", "90\\%-quantile", "Variance")) + tt_spacer_row(-5) +
    tt_rule_mid() +
    tt_numeric_row(100*as.numeric(res_means[c("mean","median","d1","d9","var")]),dec=2) +
    tt_spacer_row(-9) +
    tt_numeric_row(100*as.numeric(res_means_sd[c("mean","median","d1","d9","var")]),dec=2,surround="{\\scriptsize (%s)}") +
    tt_rule_bottom()

  tab = tt_tabularize(tt,header = "ccccc", pretty_rules=F)
  tt_save(tab,filename="inst/figures/tab-dynamic-main.tex",stand_alone=F)
  tt_save(tab,filename="~/tmp/tmp.tex",stand_alone=T)


  # SECOND TABLE
  tt =
    tt_rule_top() +
    TR("") %:% tt_text_row("{ \\bf Variance decomposition ($\\times 100$) }",cspan = 5,center = "c") + tt_spacer_row(-3) +
    tt_text_row(c("","$\\frac{Var(\\alpha)}{Var(y)}$",
                  "$\\frac{Var(\\psi)}{Var(y)}$",
                  "$\\frac{2 Cov(\\alpha,\\psi)}{Var(y)}$",
                  "$\\frac{Var(\\varepsilon)}{Var(y)}$",
                  "$Corr(\\alpha,\\psi)$")) +
    tt_rule_mid() +
    TR("estimate") %:% tt_numeric_row(100*as.numeric(main[vnames]),percent=F,dec=2) +
    tt_rule_mid_partial(list(c(2,6))) +
    TR("bootstrap mean")             %:% TR(100*as.numeric(bs_vdec_m1[vnames]),percent=F,dec=2)  +
    TR("bootstrap 2.5\\%-quantile")  %:% TR(100*as.numeric(bs_vdec_q0[vnames]),percent=F,dec=2)  +
    TR("bootstrap 97.5\\%-quantile") %:% TR(100*as.numeric(bs_vdec_q1[vnames]),percent=F,dec=2)  +
    tt_spacer_row(6) +
    TR("") %:% tt_text_row(" { \\bf Reallocation exercise ($\\times 100$) }",cspan = 5,center = "c") + tt_spacer_row(-5) +
    TR("") %:% tt_text_row( c("Mean","Median","10\\%-quantile", "90\\%-quantile", "Variance")) + tt_spacer_row(-5) +
    tt_rule_mid() +
    TR("estimate") %:% TR(100*as.numeric(res_means[c("mean","median","d1","d9","var")]),dec=2)  +
    tt_rule_mid_partial(list(c(2,6))) +
    TR("bootstrap mean")             %:% TR(100*as.numeric(res_means_m1[c("mean","median","d1","d9","var")]),dec=2)  +
    TR("bootstrap 2.5\\%-quantile")  %:% TR(100*as.numeric(res_means_q0[c("mean","median","d1","d9","var")]),dec=2)  +
    TR("bootstrap 97.5\\%-quantile") %:% TR(100*as.numeric(res_means_q1[c("mean","median","d1","d9","var")]),dec=2)  +
    tt_rule_bottom()

  tab = tt_tabularize(tt,header = "lccccc", pretty_rules=F)
  tt_save(tab,filename="inst/figures/tab-dynamic-main-bootstrap-details.tex",stand_alone=F)

}



#' plotting the means and proportions for the main estimate
#' of the mixture model.
fig.static.mixt.means.residuals <- function() {
  res_main      = res.load("m2-mixt-d2003-main-residuals")
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  I  = order(colSums(-res_main$model$A1))
  If = order(rowSums(res_main$model$A1*res_main$model$pk0[1,,]))

  gp = wplot(res_main$model$A1[If,I]) + xlab("firm class k") +theme_bw() +
    scale_y_continuous("log-earnings")+theme(legend.position = "none") +
    scale_color_brewer(palette="Spectral")
  ggsave("inst/figures/fig-mixt2-d2003-resid-wage.pdf",gp,width = 6.5,height = 5)

  # proportions
  dpk1 = m2.get.pk1(res_main$model)
  pk_m  = acast(dpk1[,sum(pr_j1j2k),list(j1,k)],j1~k,value.var = "V1")
  NNs = res_main$model$NNs*10 # used 10% sample
  NNm = res_main$model$NNm
  share_s  = sum(NNs)/(sum(NNm) + sum(NNs))
  pk_unc = share_s*rdim(res_main$model$pk0[,,I],res_main$model$nf,res_main$model$nk) +
    (1- share_s) * pk_m

  dd = data.table(melt(pk_unc,c('l','k')))
  Ir = rank(rowSums(res_main$model$A1*res_main$model$pk0[1,,]))
  dd[,l := Ir[l]]
  gp = ggplot(dd,aes(x=factor(l),y=value,fill=factor(k))) +
    geom_bar(position="stack",stat = "identity") + theme_bw() +
    theme(legend.position = "none") +
    xlab("firm class k") + ylab("type proportions") +
    scale_fill_brewer(palette="Spectral")

  ggsave("inst/figures/fig-mixt2-d2003-resid-pk.pdf",gp,width = 6.5,height = 5)
}


fig.static.mixt.means.withx <- function() {
  res_wx = res.load("m2-mixt-d2003-main-withx")
  I = order(colSums(-res_wx$mix$model$A1))

  dd = data.table(melt(res_wx$mix$model$pk0[,,I],c('x','l','k')))
  dd[,educ:= sprintf("educ %i",((x-1) %%3)+1)]
  dd[,age := sprintf("age %i",ceiling(x/3)) ]
  gp = ggplot(dd,aes(x=factor(l),y=value,fill=factor(k))) + geom_bar(position="stack",stat = "identity") + theme_bw() +
    facet_grid(educ~age) +
    theme(legend.position = "none") +
    xlab("firm class k") + ylab("type proportions") +
    scale_fill_brewer(palette="Spectral")

  ggsave("inst/figuresfig-mixt2-d2003-pk-withx.pdf",gp,width = 6.5,height = 5)
}

#' plotting the means and proportions for the split sample
#' estimates
fig.static.mixt.splits <- function() {
  res_split = res.load("m2-mixt-y2003-splits")

  I = order(-colSums(res_split$split1$model$A1))
  gp = wplot(res_split$split1$model$A1[,I]) +  xlab("firm class k") +theme_bw() +
    scale_y_continuous("log-earnings")+theme(legend.position = "none") +
    scale_color_brewer(palette="Spectral")
  ggsave("inst/figuresfig-mixt2-d2003-S1-wage.pdf",gp,width = 6.5,height = 5)

  dpk1 = m2.get.pk1(res_split$split1$model)
  pk_m  = acast(dpk1[,sum(pr_j1j2k),list(j1,k)],j1~k,value.var = "V1")
  NNs = res_split$split1$model$NNs*10 # used 10% sample
  NNm = res_split$split1$model$NNm
  share_s  = sum(NNs)/(sum(NNm) + sum(NNs))
  pk_unc = share_s*rdim(res_split$split1$model$pk0[,,I],res_split$split1$model$nf,res_split$split1$model$nk) +
    (1- share_s) * pk_m

  gp = pplot(pk_unc)+
    theme(legend.position = "none") +  xlab("firm class k") + ylab("type proportions") +
    scale_fill_brewer(palette="Spectral")
  ggsave("inst/figuresfig-mixt2-d2003-S1-pk.pdf",gp,width = 6.5,height = 5)

  I = order(-colSums(res_split$split2$model$A1))
  gp = wplot(res_split$split2$model$A1[,I]) +  xlab("firm class k") +theme_bw() +
    scale_y_continuous("log-earnings")+theme(legend.position = "none")  +
    scale_color_brewer(palette="Spectral")
  ggsave("inst/figuresfig-mixt2-d2003-S2-wage.pdf",gp,width = 6.5,height = 5)

  dpk1 = m2.get.pk1(res_split$split2$model)
  pk_m  = acast(dpk1[,sum(pr_j1j2k),list(j1,k)],j1~k,value.var = "V1")
  NNs = res_split$split2$model$NNs*10 # used 10% sample
  NNm = res_split$split2$model$NNm
  share_s  = sum(NNs)/(sum(NNm) + sum(NNs))
  pk_unc = share_s*rdim(res_split$split2$model$pk0[,,I],res_split$split2$model$nf,res_split$split2$model$nk) +
    (1- share_s) * pk_m

  gp = pplot(pk_unc)+
    theme(legend.position = "none") + xlab("firm class k") + ylab("type proportions") +
    scale_fill_brewer(palette="Spectral")
  ggsave("inst/figuresfig-mixt2-d2003-S2-pk.pdf",gp,width = 6.5,height = 5)
}

#' this plots the connectedness versus likelood for the 50
#' starting values in the main estimatation.
fig.static.mixt.connectedness <- function() {

  res_main      = res.load("m2-mixt-y2003-main-fixb")
  res_reps = res_main$second_stage_reps

  gp = ggplot(res_reps,aes(x=lik_mixt,y=connectedness,shape=factor(sel),color=factor(sel))) +
    geom_point(size=5) +theme_bw() + geom_point(data=res_reps[sel==1],size=10) +
    xlab("likelihood") +
    scale_shape_manual(values=c(20,17,8)) +
      theme(legend.position = "none")

  ggsave("inst/figuresfig-mixt2-d2003-startings.pdf",gp,width = 6.5,height = 5)
}

#' ploits the earnings before and after a move
fig.static.y1plusy2_k1k2 <- function() {
  mstats        = res.load("m2-stats-d2003")$mstats

  rr1 = mstats[,list(m12=(m1+m2)/2,N), list(j1,j2)]
  rr2 = mstats[,list(m12=(m1+m2)/2), list(j1=j2,j2=j1)]
  setkey(rr1,j1,j2)
  setkey(rr2,j1,j2)
  rr1[, m12b := rr2[rr1,m12]]

  ggplot(rr1[j1<j2],aes(x=m12,y=m12b,size=N)) + geom_point() +
    theme_bw() + geom_abline(linetype=2) +
    xlab("mean log earnings, upward mobility") +
    ylab("mean log earnings, downard mobility") +
    theme(legend.position = "none")

  # computes share of movers moving up
  rr1[j1<j2][m12b>m12,sum(N)]/rr1[j1<j2][,sum(N)]

  ggsave("inst/figuresfig-ymobility.pdf",dpi=100,width = 8,height = 7)
}

fig.static.mixture.fit <- function() {
  m2.fit = res.load("m2-mixt-d2003-main-fit")

  gp = ggplot(m2.fit$dd_m,aes(x=x,y=y1)) + geom_line() + geom_line(aes(y=y2),color="red",linetype=2) +
    facet_grid(j1b~j2b) + theme_bw()  +
    scale_x_continuous("log earnings",breaks=c(9.0,10.0,11.0)) + ylab("period 1 density")
  ggsave("inst/figuresfig-m2-fit-density-movers.pdf",gp,dpi=100,width = 6,height = 5)

  gp = ggplot(m2.fit$dd_s,aes(x=x,y=y1)) + geom_line() + geom_line(aes(y=y2),color="red",linetype=2) +
    facet_wrap(~j1,nrow=2) + theme_bw()  +
    xlab("log earnings") + ylab("period 1 density")
  ggsave("inst/figuresfig-m2-fit-density-stayers.pdf",gp,dpi=100,width = 7.5,height = 5)

  # fit ov mean/covariance
  dd1 = data.table(m2.fit$rr_s)
  dd2 = data.table(m2.fit$rr_m)
  dd1[,variable:=sprintf("stayer-%s",variable)]
  dd1[,variable:=sprintf("movers-%s",variable)]
  ggplot(dd[variable %in% c("stayers-m1","stayers-v1","movers-cov12")],aes(x=data,y=imp)) + geom_point() +
    facet_wrap(~variable,scale="free") + geom_abline() + theme_bw()

  res_main = res.load("m2-mixt-y2003-main-fixb")
  NNm = res_main$model$NNm

  dd = data.table(m2.fit$rr_m)
  dd[,N := NNm[j1,j2],list(j1,j2)]
  gp = ggplot(dd[variable %in% c("cov12")],aes(x=data,y=imp,size=N)) + geom_point(alpha=0.4) +
    geom_abline() + theme_bw() + ylim(-0.2,0.3) + xlab("covariance data") + ylab("covariance model") +
    theme(legend.position = "none")

  ggsave("inst/figuresfig-m2-fit-density-cov.pdf",gp,dpi=100,width = 5,height = 5)

}

fig.static.mixture.mobility <- function() {
  #res_main      = res.load("m2-mixt-y2003-main-fixb")
  res_main   = res.load("m4-mixt-d2003-main")
  I = rank(colSums(res_main$model$A1))

  dpk1 = m2.get.pk1(res_main$model)
  dpk1[,k := I[k]]
  dpk1[, pr_tr := pr_j1j2k/sum(pr_j1j2k),list(k)]
  dpk1[,wn :=  sprintf("worker type %i",k) ]

  gp = ggplot(dpk1,aes(x=factor(j1))) +
    geom_abline(linetype=2,alpha=0.3) +
    geom_point(aes(y=factor(j2),alpha=pr_tr,size=pr_tr)) + facet_wrap(~wn,nrow=2) + theme_bw() +
    xlab("firm class k, period 1") + ylab("firm class k', period 2") +
    theme(legend.position = "none")
  ggsave("inst/figuresfig-mixture-dynamic-mobility-joint.pdf",gp,dpi=100,width = 7.5,height = 5)

  #dpk1[, pr_tr := pr_j1j2k/sum(pr_j1j2k),list(j1,k)]
  dpkb = dpk1[,wtd.mean(j2,pr_tr),list(j1,k)]
  ggplot(dpk1,aes(x=factor(j1))) +
    geom_abline(linetype=2,alpha=0.3) +
    facet_wrap(~k,nrow=2) + theme_bw() +
    geom_step(aes(y=V1,group=k),size=1,alpha=1,data=dpkb) +
    scale_y_continuous("mean firm class k', period 2",breaks=1:10,minor_breaks = 1:10) +
    coord_cartesian(xlim=c(1,10),ylim=c(1,10)) +
    xlab("firm class k, period 1")


  dpkb = dpk1[,sum(pr_tr[j2>j1])/sum(pr_tr),list(j1,k)]
  ggplot(dpkb,aes(x=factor(j1))) +
     theme_bw() +
    geom_line(aes(y=V1,group=k)) +
    scale_y_continuous("mean firm class k', period 2",breaks=1:10,minor_breaks = 1:10) +
    coord_cartesian(xlim=c(1,10),ylim=c(1,10)) +
    xlab("firm class k, period 1")


}


# ========== AKM ============

#' AKM plot with bias correction
fig.akm <- function() {
  res_akm  = res.load("akm-fe-bc-d2003")

  gp = ggplot(res_akm,aes(x=nm)) + geom_line(aes(y=v0)) + geom_line(aes(y=vh)) +
    geom_point(aes(y=v0),shape=1,size=3) + geom_point(aes(y=vh),shape=2,size=3) + theme_bw() +
    geom_hline(yintercept = 0.14*(0.0342*0.75) ,linetype=2) + ylim(0,0.012)  +
    ylab("variance of firm effect") + xlab("minimum number of movers per firm") +
    scale_x_continuous(breaks=c(5,10,15,20,25,30,35,40),limits = c(4,40))

  ggsave("inst/figuresfig-akm-fe.pdf",gp,dpi=100,width = 6,height = 5)
}

# =========================  DYNAMIC -- MIXTURE =====================

#' this plots the connectedness versus likelood for the 50
#' starting values in the main estimatation.
fig.dynamic.mixt.connectedness <- function() {

  res_main   = res.load("m4-mixt-d2003-main")
  res_reps = res_main$second_stage_reps

  gp = ggplot(res_reps,aes(x=lik_mixt,y=connectedness,shape=factor(sel),color=factor(sel))) +
    geom_point(size=5) +theme_bw() + geom_point(data=res_reps[sel==1],size=10) +
    xlab("likelihood") +
    scale_shape_manual(values=c(20,17,8)) +
    theme(legend.position = "none")

  ggsave("inst/figuresfig-mixt4-d2003-startings.pdf",gp,width = 6.5,height = 5)
}

#' plotting the means and proportions for the main estimate
#' of the mixture model.
fig.dynamic.mixt.means <- function() {
  res_main   = res.load("m4-mixt-d2003-main")
  model_mixt_bs  = m4.getboostrap()$mixt_all
  flog.info("using %i bootstrap replications", length(model_mixt_bs))

  cstats        = res.load("m4-stats-d2003")$cstats
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  res_mixt_bs = data.table(ldply(model_mixt_bs,function(x) {
    I = order(colSums(-x$model$A2s))
    #If = order(rowSums(x$model$A1*x$model$pk0[,]))
    melt(x$model$A2s[,I],c('l','k')) }
  ))
  res_sd      = res_mixt_bs[,list(q0=quantile(value,0.025),q1=quantile(value,0.975),m=mean(value)),list(l,k)]

  I = order(colSums(-res_main$model$A2s))
  res_sd[,value:=res_main$model$A2s[l,I[k]] ,list(k,l)]

  lmin = min(cstats[,m1-3*sd1])
  lmax = max(cstats[,m1+3*sd1])
  cstats$k=1

  wplot(res_main$model$A2s[,I]) + xlab("firm class k") +theme_bw() +
    scale_y_continuous("log-earnings")+theme(legend.position = "none") +
    geom_errorbar(aes(ymin=q0,ymax=q1),data=res_sd,width=0.2) +
    coord_cartesian(ylim = c(lmin,lmax)) #+
  #geom_ribbon(aes(ymin=m1-2*sd1,ymax=m1+2*sd1,x=j1,y=m1),data=cstats,alpha=0.1,color="grey")

  gp = wplot(res_main$model$A2s[,I]) + xlab("firm class k") +theme_bw() +
    scale_y_continuous("log-earnings")+theme(legend.position = "none") +
    geom_errorbar(aes(ymin=value-(q1-q0)/2,ymax=value+(q1-q0)/2),data=res_sd,width=0.2) +
    #coord_cartesian(ylim = c(lmin,lmax)) +
    #geom_ribbon(aes(ymin=m1-2*sd1,ymax=m1+2*sd1,x=j1,y=m1),data=cstats,alpha=0.05,color=NA)+
    scale_color_brewer(palette="Spectral")

  ggsave("inst/figures/fig-mixt4-d2003-wage.pdf",gp,width = 6.5,height = 5)
  flog.info("creating fig-mixt4-d2003-wage.png")

  gp = wplot(res_main$model$A2s[,I]) + xlab("firm class k") +theme_bw() +
    scale_y_continuous("log-earnings")+theme(legend.position = "none") +
    geom_errorbar(aes(ymin=q0,ymax=q1),data=res_sd,width=0.2) +
    #coord_cartesian(ylim = c(lmin,lmax)) +
    # geom_ribbon(aes(ymin=m1-2*sd1,ymax=m1+2*sd1,x=j1,y=m1),data=cstats,alpha=0.05,color=NA) +
    scale_color_brewer(palette="Spectral")

  ggsave("inst/figures/fig-mixt4-d2003-wage-uc.pdf",gp,width = 6.5,height = 5)
  flog.info("creating fig-mixt4-d2003-wage-uc.png")

  # proportions
  dpk1 = m2.get.pk1(res_main$model)
  pk_m  = acast(dpk1[,sum(pr_j1j2k),list(j1,k)],j1~k,value.var = "V1")
  NNs = res_main$model$NNs*10 # used 10% sample
  NNm = res_main$model$NNm
  share_s  = sum(NNs)/(sum(NNm) + sum(NNs))
  pk_unc = share_s*rdim(res_main$model$pk0[,I],res_main$model$nf,res_main$model$nk) +
    (1- share_s) * pk_m

  dd = melt(pk_unc,c('l','k'))
  gp = ggplot(dd,aes(x=factor(l),y=value,fill=factor(k))) +
    geom_bar(position="stack",stat = "identity") + theme_bw() +
    theme(legend.position = "none") +
    xlab("firm class k") + ylab("type proportions") +
    scale_fill_brewer(palette="Spectral")

  ggsave("inst/figuresfig-mixt4-d2003-pk.pdf",gp,width = 6.5,height = 5)
  flog.info("creating fig-mixt4-d2003-pk.png")
}


#' plotting the means and proportions for the split sample
#' estimates
fig.dynamic.mixt.splits <- function() {
  res_split  = res.load("m4-mixt-d2003-splits")

  I = order(-colSums(res_split$split1$mixt$model$A2s))
  gp = wplot(res_split$split1$mixt$model$A2s[,I]) +  xlab("firm class k") +theme_bw() +
    scale_y_continuous("log-earnings")+theme(legend.position = "none") +
    scale_color_brewer(palette="Spectral")
  ggsave("inst/figuresfig-mixt4-d2003-S1-wage.pdf",gp,width = 6.5,height = 5)

  dpk1 = m2.get.pk1(res_split$split1$mixt$model)
  pk_m  = acast(dpk1[,sum(pr_j1j2k),list(j1,k)],j1~k,value.var = "V1")
  NNs = res_split$split1$mixt$model$NNs*10 # used 10% sample
  NNm = res_split$split1$mixt$model$NNm
  share_s  = sum(NNs)/(sum(NNm) + sum(NNs))
  pk_unc = share_s*rdim(res_split$split1$mixt$model$pk0[,I],res_split$split1$mixt$model$nf,res_split$split1$mixt$model$nk) +
    (1- share_s) * pk_m

  gp = pplot(pk_unc)+
    theme(legend.position = "none") +  xlab("firm class k") + ylab("type proportions") +
    scale_fill_brewer(palette="Spectral")
  ggsave("inst/figuresfig-mixt4-d2003-S1-pk.pdf",gp,width = 6.5,height = 5)

  # I = order(-colSums(res_split$split2$model$A1))
  # gp = wplot(res_split$split2$model$A1[,I]) +  xlab("firm class k") +theme_bw() +
  #   scale_y_continuous("log-earnings")+theme(legend.position = "none") +
  #   scale_color_brewer(palette="Spectral")
  # ggsave("inst/figuresfig-mixt2-d2003-S2-wage.png",gp,width = 6.5,height = 5)
  #
  # dpk1 = m2.get.pk1(res_split$split2$model)
  # pk_m  = acast(dpk1[,sum(pr_j1j2k),list(j1,k)],j1~k,value.var = "V1")
  # NNs = res_split$split2$model$NNs*10 # used 10% sample
  # NNm = res_split$split2$model$NNm
  # share_s  = sum(NNs)/(sum(NNm) + sum(NNs))
  # pk_unc = share_s*rdim(res_split$split2$model$pk0[,,I],res_split$split2$model$nf,res_split$split2$model$nk) +
  #   (1- share_s) * pk_m
  #
  # gp = pplot(pk_unc)+
  #   theme(legend.position = "none") +
  #   scale_fill_brewer(palette="Spectral")
  # ggsave("inst/figuresfig-mixt2-d2003-S2-pk.png",gp,width = 6.5,height = 5)
}


fig.dynamic.mixture.fit <- function() {
  m2.fit = res.load("m4-mixt-d2003-main-fit")

  gp = ggplot(m2.fit$dd_m,aes(x=x,y=y1)) + geom_line() + geom_line(aes(y=y2),color="red",linetype=2) +
    facet_grid(j1b~j2b) + theme_bw()  +
    scale_x_continuous("log earnings",breaks=c(9.0,10.0,11.0)) + ylab("period 2 density")
  ggsave("inst/figuresfig-m4-fit-density-movers.pdf",gp,dpi=100,width = 6,height = 5)

  gp = ggplot(m2.fit$dd_s,aes(x=x,y=y1)) + geom_line() + geom_line(aes(y=y2),color="red",linetype=2) +
    facet_wrap(~j1,nrow = 2) + theme_bw()  +
    xlab("log earnings") + ylab("period 2 density")
  ggsave("inst/figuresfig-m4-fit-density-stayers.pdf",gp,dpi=100,width = 6,height = 5)

  # fit ov mean/covariance
  dd1 = data.table(m2.fit$rr_s)
  dd2 = data.table(m2.fit$rr_m)
  dd1[,variable:=sprintf("stayer-%s",variable)]
  dd1[,variable:=sprintf("movers-%s",variable)]
  ggplot(dd[variable %in% c("stayers-m1","stayers-v1","movers-cov12")],aes(x=data,y=imp)) + geom_point() +
    facet_wrap(~variable,scale="free") + geom_abline() + theme_bw()

  res_main   = res.load("m4-mixt-d2003-main")
  NNm = res_main$model$NNm

  dd = data.table(m2.fit$rr_m)
  dd[,N := NNm[j1,j2],list(j1,j2)]
  gp = ggplot(dd[variable %in% c("cov12")],aes(x=data,y=imp,size=N)) + geom_point(alpha=0.4) +
    geom_abline() + theme_bw() + ylim(-0.2,0.3) + xlab("covariance data") + ylab("covariance model") +
    theme(legend.position = "none")

  ggsave("inst/figuresfig-m4-fit-density-cov.pdf",gp,dpi=100,width = 5,height = 5)

  gp = ggplot(m2.fit$dd_m_g,aes(x=x,y=y1)) + geom_line() + geom_line(aes(y=y2),color="red",linetype=2) +
    facet_grid(j1b~j2b) + theme_bw()  +
    scale_x_continuous("log earnings",breaks=c(-0.5,0.0,0.5)) + ylab("density of wage growth (t=1,2), stayers")
  ggsave("inst/figures/fig-m4-fit-density-movers-growth.pdf",gp,dpi=100,width = 6,height = 5)

  gp = ggplot(m2.fit$dd_s_g,aes(x=x,y=y1)) + geom_line() + geom_line(aes(y=y2),color="red",linetype=2) +
    facet_wrap(~j1,nrow = 2) + theme_bw()  +
    xlab("log earnings") + ylab("earning growth density")
  ggsave("inst/figures/fig-m4-fit-density-stayers-growth.pdf",gp,dpi=100,width = 6,height = 5)



}

# plots the slice of the rho
# from the covariance structure
fig.dyn.rho <-function() {
  path_archive= "inst/fig/res/"
  rho_evals = data.table(archive.get("rhos_grid","res-2003-dynamic.dat"))
  rho_evals[name=="rho1",name:="rho 1|2"]
  rho_evals[name=="rho4",name:="rho 4|3"]
  ggplot(rho_evals[name!="rho32"],aes(x=x,y=y)) + geom_line() + facet_grid(~name) +
    theme_bw() + geom_vline(aes(xintercept=best),linetype=2) +
    coord_cartesian(ylim = c(0,40)) + xlab("value of rho") + ylab("value of the objective")
  ggsave("inst/figuresfig-rho-slice.pdf",dpi=100,width = 8,height = 4)
}

fig.dyn.rho.fit <-function() {
  res_main   = res.load("m4-mixt-d2003-main")
  res_rhos = res_main$res_rhos
  ddres = data.table(data=res_rhos$dep,model=res_rhos$fitted,w=res_rhos$weights)

  ggplot(ddres,aes(x=data,y=model,size=w^4)) +geom_point(alpha=0.3) +
    geom_abline(linetype=2)  + theme_bw() +
    theme(legend.position = "none")
  ggsave("inst/figuresfig-rho-fit.pdf",dpi=100,width = 5,height = 4)

  res_main_rho   = res.load("m4-mixt-d2003-rho-check")
  res_rhos = res_main_rho$res_rhos
  ddres = data.table(data=res_rhos$dep,model=res_rhos$fitted,w=res_rhos$weights)

  ggplot(ddres,aes(x=data,y=model,size=w^4)) +geom_point(alpha=0.3) +
    geom_abline(linetype=2)  + theme_bw() +
    theme(legend.position = "none")
  ggsave("inst/figuresfig-rho-diff-fit.pdf",dpi=100,width = 5,height = 4)
}

table.m4.parameters <- function() {

  res_main   = res.load("m4-mixt-d2003-main")
  res_bs     = m4.getboostrap()

  model_mixt_bs = data.table(ldply(res_bs$mixt_all,function(x) {
    I = order(colSums(-x$model$A2s))
    #If = order(rowSums(x$model$A2s*x$model$pk0[,]))
    data.frame(k=1:10,A2mb=x$model$A2mb,A3mb=x$model$A3mb) }
  ))

  model_mixt_bs2 = data.table(ldply(res_bs$mixt_all,function(x) {
    data.frame(x$model[c("B32m","B32s","B12","B43")]) }
  ))

  model_mixt_bs  = model_mixt_bs[,list(A2mb = sd(A2mb),A3mb = sd(A3mb)),list(k)]
  model_mixt_bs2 = model_mixt_bs2[,list(B32m = sd(B32m),B32s = sd(B32s),B12 = sd(B12),B43 = sd(B43) )]

  tt =
    tt_rule_top() +
    TR(c("","Earnings effects $\\xi_2(k')$ of future firm classes"),cspan=c(1,9)) +
    tt_rule_mid_partial(list(c(2,10))) +
    TR(c("$k'=$",paste(2:10))) +
    TR("estimate") %:% TR(res_main$model$A2mb[2:10],dec=3)  + tt_spacer_row(-8) +
    TR("")         %:% TR(model_mixt_bs$A2mb[2:10],dec=3,surround="{\\footnotesize (%s)}") +
    tt_rule_mid() +
    TR(c("","Earnings effects $\\xi_3(k)$ of past firm classes"),cspan=c(1,9)) +
    tt_rule_mid_partial(list(c(2,10))) +
    TR(c("$k=$",paste(2:10))) +
    TR("estimate") %:% TR(res_main$model$A3mb[2:10],dec=3)  + tt_spacer_row(-8) +
    TR("")         %:% TR(model_mixt_bs$A3mb[2:10],dec=3,surround="{\\footnotesize (%s)}") +
    tt_rule_mid() +
    TR(c("","Persistence parameters $\\rho$"),cspan=c(1,9)) +
    tt_rule_mid_partial(list(c(2,10))) +
    TR(c("","","$\\rho_{1|2}$","","$\\rho_{3|2}^m$","","$\\rho_{3|2}^s$","","$\\rho_{4|3}$")) +
    TR(" ") %:%  TR(" ") %:% TR(res_main$model$B12,dec=3) %:% TR("") %:%
                     TR(res_main$model$B32m,dec=3) %:% TR("") %:%
                     TR(res_main$model$B32s,dec=3) %:% TR("") %:%
                     TR(res_main$model$B43,dec=3)  + tt_spacer_row(-8) +
    TR(" ") %:%  TR(" ") %:% TR(model_mixt_bs2$B12,dec=3,surround="{\\footnotesize (%s)}") %:% TR("") %:%
    TR(model_mixt_bs2$B32m,dec=3,surround="{\\footnotesize (%s)}") %:% TR("") %:%
    TR(model_mixt_bs2$B32s,dec=3,surround="{\\footnotesize (%s)}") %:% TR("") %:%
    TR(model_mixt_bs2$B43,dec=3,surround="{\\footnotesize (%s)}") +
    tt_rule_bottom()

    tab = tt_tabularize(tt,"l ccccccccc")
    tt_save(tab,filename="inst/figures/tab-dyn-params.tex",stand_alone=F)


}



table.endogeneousMobility <- function() {
  rr = res.load("m4-cf-mobility")

  tt =
    tt_rule_top() +
    TR(c(""," \\bf All"),cspan=c(1,4)) +
    TR(c("","All movers","$k'{=}1{-}3$", "$k'{=}4{-}7$", "$k'{=}8{-}10$")) +
    tt_rule_mid_partial(list(c(2,2),c(3,5))) +
    TR("$k{=}1{-}3$")  %:% TR(100*rr$main[j1=="k=1..3"][cond==2,value],dec=2) + tt_spacer_row(-8) +
    TR("")         %:% TR(100*rr$stats[j1=="k=1..3"] [cond==2,sd],dec=2,surround="{\\footnotesize (%s)}") +
    TR("$k{=}4{-}7$")  %:% TR(100*rr$main[j1=="k=4..7"]  [cond==2,value],dec=2) + tt_spacer_row(-8) +
    TR("")         %:% TR(100*rr$stats[j1=="k=4..7"] [cond==2,sd],dec=2,surround="{\\footnotesize (%s)}") +
    TR("$k{=}8{-}10$") %:% TR(100*rr$main[j1=="k=8..10"] [cond==2,value],dec=2) + tt_spacer_row(-8) +
    TR("")         %:% TR(100*rr$stats[j1=="k=8..10"][cond==2,sd],dec=2,surround="{\\footnotesize (%s)}") +
    tt_rule_mid() +
    TR(c(""," \\bf First conditional decile of earnings"),cspan=c(1,4)) +
    TR(c("","All movers","$k'{=}1{-}3$", "$k'{=}4{-}7$", "$k'{=}8{-}10$")) +
    tt_rule_mid_partial(list(c(2,2),c(3,5))) +
    TR("$k{=}1{-}3$")  %:% TR(100*rr$main[j1=="k=1..3"]  [cond==1,value],dec=2) + tt_spacer_row(-8) +
    TR("")         %:% TR(100*rr$stats[j1=="k=1..3"] [cond==1,sd],dec=2,surround="{\\footnotesize (%s)}") +
    TR("$k{=}4{-}7$")  %:% TR(100*rr$main[j1=="k=4..7"]  [cond==1,value],dec=2) + tt_spacer_row(-8) +
    TR("")         %:% TR(100*rr$stats[j1=="k=4..7"] [cond==1,sd],dec=2,surround="{\\footnotesize (%s)}") +
    TR("$k{=}8{-}10$") %:% TR(100*rr$main[j1=="k=8..10"] [cond==1,value],dec=2) + tt_spacer_row(-8) +
    TR("")         %:% TR(100*rr$stats[j1=="k=8..10"][cond==1,sd],dec=2,surround="{\\footnotesize (%s)}") +
    tt_rule_mid() +
    TR(c(""," \\bf Tenth conditional decile of earnings"),cspan=c(1,4)) +
    TR(c("","All movers","$k'{=}1{-}3$", "$k'{=}4{-}7$", "$k'{=}8{-}10$")) +
    tt_rule_mid_partial(list(c(2,2),c(3,5))) +
    TR("$k{=}1{-}3$")  %:% TR(100*rr$main[j1=="k=1..3"]  [cond==3,value],dec=2) + tt_spacer_row(-8) +
    TR("")         %:% TR(100*rr$stats[j1=="k=1..3"] [cond==3,sd],dec=2,surround="{\\footnotesize (%s)}") +
    TR("$k{=}4{-}7$")  %:% TR(100*rr$main[j1=="k=4..7"]  [cond==3,value],dec=2) + tt_spacer_row(-8) +
    TR("")         %:% TR(100*rr$stats[j1=="k=4..7"] [cond==3,sd],dec=2,surround="{\\footnotesize (%s)}") +
    TR("$k{=}8{-}10$") %:% TR(100*rr$main[j1=="k=8..10"] [cond==3,value],dec=2) + tt_spacer_row(-8) +
    TR("")         %:% TR(100*rr$stats[j1=="k=8..10"][cond==3,sd],dec=2,surround="{\\footnotesize (%s)}") +
    tt_rule_top()

  tab = tt_tabularize(tt,header = "l cccc", pretty_rules=F)
  tt_save(tab,filename="inst/figures/tab-dyn-endogenous-mobility.tex",stand_alone=F)
}

table.statedependence <- function() {

  rr = res.load("m4-cf-state-dependence")
  r1 =as.list(rr$main)
  names(r1) = rr$variable
  r2 = as.list(rr$sd)
  names(r2) = rr$variable

  tt =
    tt_rule_top() +
    TR(c("","","within current firm variance"),cspan=c(1,1,2)) +
    tt_rule_mid_partial(list(c(3,4))) +
    TR(c("","between current firm","between past firm","residual")) +
    tt_rule_mid_partial(list(c(2,2),c(3,4))) +
    TR("year after the move (2004)")     %:% TR(100*c(r1$y3_Wa_Bk2,r1$y3_Wak2_Bk1,r1$y3_Wak1k2),dec=2)+tt_spacer_row(-8) +
    TR("")                               %:% TR(100*c(r2$y3_Wa_Bk2,r2$y3_Wak2_Bk1,r2$y3_Wak1k2),dec=2,surround="{\\footnotesize (%s)}") +
    TR("two year after the move (2005)") %:% TR(100*c(r1$y4_Wa_Bk2,r1$y4_Wak2_Bk1,r1$y4_Wak1k2),dec=2)+tt_spacer_row(-8) +
    TR("")                               %:% TR(100*c(r2$y4_Wa_Bk2,r2$y4_Wak2_Bk1,r2$y4_Wak1k2),dec=2,surround="{\\footnotesize (%s)}") +
    tt_rule_mid() +
    TR(c("","between past firm variance"),cspan=c(1,3)) +
    tt_rule_mid_partial(list(c(2,4))) +
    TR(c("","total","network effect","state dependence")) +
    tt_rule_mid_partial(list(c(2,4))) +
    TR("year after the move (2004)")     %:% TR(100*c(r1$y3_Wa_Bk1,r1$y3net_Wa_Bk1,r1$y3var_diff),dec=2)+tt_spacer_row(-8) +
    TR("")                               %:% TR(100*c(r2$y3_Wa_Bk1,r2$y3net_Wa_Bk1,r2$y3var_diff),dec=2,surround="{\\footnotesize (%s)}") +
    TR("two year after the move (2005)") %:% TR(100*c(r1$y4_Wa_Bk1,r1$y4net_Wa_Bk1,r1$y4var_diff),dec=2)+tt_spacer_row(-8) +
    TR("")                               %:% TR(100*c(r2$y4_Wa_Bk1,r2$y4net_Wa_Bk1,r2$y4var_diff),dec=2,surround="{\\footnotesize (%s)}") +
    tt_rule_bottom()

  tab = tt_tabularize(tt,header = "l ccc", pretty_rules=F)
  tt_save(tab,filename="inst/figures/tab-dyn-state-dependence.tex",stand_alone=F)



}


# ============================ SHIMER SMITH ==================


fig.shimersmith.model <- function() {
  r = res.load("m2-shimersmith-mc")
  # make figure
  model.pam = r$pam_6x10$model
  model.nam = r$nam_6x10$model
  plotF1  <- wireframe(value ~ x* y, melt(model.pam$F,c('x','y')),zlab="",main="Production PAM"  )
  plotH1  <- wireframe(value ~ x* y, melt(model.pam$H,c('x','y')),zlab="",main="Allocation PAM"    )
  plotS1  <- wireframe(value ~ x* y, melt(model.pam$S,c('x','y')),zlab="",main="Surplus PAM"  )
  plotM1  <- levelplot(value ~ x* y, melt(ifelse(model.pam$S<=0,-Inf,model.nam$S),c('x','y')),zlab="",main="Surplus"  )
  plotF2  <- wireframe(value ~ x* y, melt(model.nam$F,c('x','y')),zlab="",main="Production NAM"  )
  plotH2  <- wireframe(value ~ x* y, melt(model.nam$H,c('x','y')),zlab="",main="Allocation NAM"    )
  plotS2  <- wireframe(value ~ x* y, melt(model.nam$S,c('x','y')),zlab="",main="Surplus NAM"  )
  plotM2  <- levelplot(value ~ x* y, melt(ifelse(model.nam$S<=0,-Inf,model.nam$S),c('x','y')),zlab="",main="Surplus"  )
  #grid.arrange(plotF1,plotF2,plotS1,plotS2,plotH1,plotH2, ncol=2)

  pdf("inst/figuresplot-shimer-model.pdf",width=12,height = 9)
  grid.arrange(plotF1,plotS1,plotH1,plotF2,plotS2,plotH2, ncol=3)
  dev.off()
}

fig.shimersmith.cardkline <- function() {

  p <- blm:::initp( b=0.3, c=0, sz=1, nx=6, ny = 10)
  p$pf = function(x,y,z=0,p)  ( 0.5*x^p$rho +  0.5*y^p$rho )^(1/p$rho) + p$ay # haggerdorn law manovski
  p$ay  = 0.7 #0.5
  p$rho = -3 # 2.5 #-1.5 # 2.5 # -1.5

  r = list()
  p$rho = -3 # 2.5 #-1.5 # 2.5 # -1.5
  p$lb1=0
  model <- shimersmith.solve.otj(p,maxiter=1000,enditer=1000)
  model$wage = model$wage +1  # shift everything up to have well defined logs

  # set sep to large value for simulation
  p$sep=1.0
  p$lb0=1.0
  cdata <- shimersmith.simulate.otj(model,p,1e6,rsq=1)

  # plot the wage by worker type
  rrw = melt(log(model$wage),c("k","l"))
  gp = ggplot(rrw,aes(x=factor(l),y=value,group=factor(k),color=factor(k))) +
    geom_line() + geom_point() +theme_bw() +  theme(legend.position = "none") +
    xlab("firm class") + ylab("log earnings")
  ggsave("inst/figuresfig-shimersmith-wages.pdf",gp,width = 4.5,height = 3)

  # construct the plot for movers between 2 periods
  jdata = cdata[mover>0]
  rr = jdata[,list(m2=mean(y1),m3=mean(y2)),list(t1,t2)]
  rr[,m1:=m2][,m4:=m3]
  rrm = melt(rr,id=c("t1","t2"))
  rrm$variable = paste(rrm$variable)
  setkey(rrm,variable)
  rrm[variable=="m1",period:=1]
  rrm[variable=="m2",period:=2]
  rrm[variable=="m3",period:=3]
  rrm[variable=="m4",period:=4]
  rrm[,trans := sprintf("%s-%s",t1,t2)]
  rrm[t1!=t2, ltype:= 2 ]
  rrm[t1==t2, ltype:= 1 ]
  rrm[,ltype := as.integer(ltype)]

  gp <- ggplot(rrm[t1%in%c(4,10)][t2%in%c(4,10)],aes(x=factor(period),y=value,color=interaction(t1,t2),group=interaction(t1,t2),linetype=factor(ltype) )) +
    geom_line(size=1.5) + theme_bw() + xlab("period") + ylab("log earnings") + theme(legend.position = "none")

  ggsave("inst/figuresfig-shimersmith-CK.pdf",gp,width = 4.5,height = 3)

}

fig.shimersmith.wages <- function() {

  res_ss = res.load("m2-shimersmith-mc")

  # PAM model
  model = res_ss$pam_6x10$model
  gp = wplot(t(log(model$wage))) + xlab("firm class k") +theme_bw() + geom_point() +
    scale_y_continuous("log-earnings")+theme(legend.position = "none")
  ggsave("inst/figuresfig-shimersmith-pam_6x10_model.pdf",gp,width = 6.5,height = 5)

  model = res_ss$pam_6x10$res_mixt$model
  gp = wplot(model$A2) + xlab("firm class k") +theme_bw() + geom_point() +
    scale_y_continuous("log-earnings")+theme(legend.position = "none")
  ggsave("inst/figuresfig-shimersmith-pam_6x10_estimate.pdf",gp,width = 6.5,height = 5)

  # NAM model
  model = res_ss$nam_6x10$model
  gp = wplot(t(log(model$wage))) + xlab("firm class k") +theme_bw() + geom_point() +
    scale_y_continuous("log-earnings")+theme(legend.position = "none")
  ggsave("inst/figuresfig-shimersmith-nam_6x10_model.pdf",gp,width = 6.5,height = 5)

  model = res_ss$nam_6x10$res_mixt$model
  I = order(-model$pk0[1,,6])
  gp = wplot(model$A2[I,]) + xlab("firm class k") +theme_bw() + geom_point() +
    scale_y_continuous("log-earnings")+theme(legend.position = "none")
  ggsave("inst/figuresfig-shimersmith-nam_6x10_estimate.pdf",gp,width = 6.5,height = 5)

  cstats = res_ss$pam_6x10$cstats
  gp = ggplot(melt(cstats,id.vars = c('j1','m','sd')), aes(x=factor(j1),y=value,group=variable)) +
    geom_point() + geom_line() + theme_bw() + xlab("firm class k") + scale_y_continuous("log-earnings") +
    geom_line(data=cstats,aes(y=m,group=NA),size=1.5)
  ggsave("inst/figuresfig-shimersmith-pam_6x10_distr.pdf",gp,width = 6.5,height = 5)

  cstats = res_ss$nam_6x10$cstats
  gp = ggplot(melt(cstats,id.vars = c('j1','m','sd')), aes(x=factor(j1),y=value,group=variable)) +
    geom_point() + geom_line() + theme_bw() + xlab("firm class k") + scale_y_continuous("log-earnings") +
    geom_line(data=cstats,aes(y=m,group=NA),size=1.5)
  ggsave("inst/figuresfig-shimersmith-nam_6x10_distr.pdf",gp,width = 6.5,height = 5)

  # re-simulate from the model to compute mstats
  model = res_ss$pam_6x10$model
  cdata <- shimersmith.simulate.otj(model,res_ss$pam_6x10$p,5e5,rsq=0.9)

  # we cluster the data
  flog.info("cluster data")
  cdata[,f1:=paste(f1)][,f2:=paste(f2)]
  setnames(cdata,"x","k")
  cdata$x=1
  sim = list(sdata=cdata[mover==0],jdata=cdata[mover>0])
  ms    = grouping.getMeasures(sim,"ecdf",Nw=20,y_var = "y1")
  grps  = grouping.classify.once(ms,k = 10,nstart = 1000,iter.max = 200,step=100)
  sim   = grouping.append(sim,grps$best_cluster)

  mstats = sim$jdata[,list(m1=mean(y1),sd1=sd(y1),
                           m2=mean(y2),sd2=sd(y2),
                           v12 = cov(y1,y2),.N),list(j1,j2)]

  rr1 = mstats[,list(m12=(m1+m2)/2,N), list(j1,j2)]
  rr2 = mstats[,list(m12=(m1+m2)/2), list(j1=j2,j2=j1)]
  setkey(rr1,j1,j2)
  setkey(rr2,j1,j2)
  rr1[, m12b := rr2[rr1,m12]]

  ggplot(rr1[j1<j2],aes(x=m12,y=m12b,size=N)) + geom_point() +
    theme_bw() + geom_abline(linetype=2) +
    xlab("mean log earnings, upward mobility") +
    ylab("mean log earnings, downard mobility") +
    theme(legend.position = "none")

  # computes share of movers moving up
  rr1[j1<j2][m12b>m12,sum(N)]/rr1[j1<j2][,sum(N)]

  ggsave("inst/figures/fig-shimersmith-pam_6x10_moverplot.pdf",gp,width = 5,height = 5)
}

fig.shimersmith.narrow <- function() {
  rr = res.load("m2-shimersmith-narrow")

  model = rr$model
  gp1 = wplot(t(log(model$wage))) + xlab("firm class k") +theme_bw() + geom_point() +
    scale_y_continuous("log-earnings")+theme(legend.position = "none")  +scale_color_grey()
  ggsave("inst/figures/fig-shimersmith-narrow_model.pdf",gp1,width = 6.5,height = 5)

  rr = data.table(melt(model$S,c("x","y")))
  gp3 = ggplot(rr[value>=0],aes(x=x,y=y,fill=value)) + geom_tile() + theme_bw() + ylab("firms") + xlab("workers") +
     theme(legend.position = "none") +scale_fill_grey()
  ggsave("inst/figures/fig-shimersmith-narrow_matchset.pdf",gp3,width = 5,height = 5)

  model = rr$res_mixt$model
  I = order(model$pk0[1,,6])
  gp2 = wplot(model$A2[I,]) + xlab("firm class k") +theme_bw() + geom_point() +
    scale_y_continuous("log-earnings")+theme(legend.position = "none")   +scale_color_grey()
  ggsave("inst/figures/fig-shimersmith-narrow_estimate.pdf",gp2,width = 6.5,height = 5)


  require(textables)

  vdec.extract <- function(vdec) {
    VV = vdec$cc
    r = data.frame(var_y=VV[1,1],var_k=VV[2,2],var_l=VV[3,3],var_e=VV[4,4],cov_kl=VV[2,3])
    r$cor_kl        = r$cov_kl/sqrt(r$var_k*r$var_l)
    r$var_l_tshare  = r$var_l/r$var_y
    r$var_k_tshare  = r$var_k/r$var_y
    r$var_e_tshare  = r$var_e/r$var_y
    r$cov_kl_tshare = r$cov_kl/r$var_y
    return(100*c(r$var_k_tshare,r$var_l_tshare,r$cov_kl_tshare,r$var_e_tshare,r$cor_kl))
  }

  tt = tt_rule_top() +
    tt_text_row(c("","$\\frac{Var(\\alpha)}{Var(y)}$",
                  "$\\frac{Var(\\psi)}{Var(y)}$",
                  "$\\frac{2 Cov(\\alpha,\\psi)}{Var(y)}$",
                  "$\\frac{Var(\\varepsilon)}{Var(y)}$",
                  "$Corr(\\alpha,\\psi)$")) +
    tt_rule_mid_partial(list(c(2,6))) +
    tt_text_row("DGP") %&% tt_numeric_row(vdec.extract(rr$vdec_true),perc=F,dec=1) +
    tt_text_row("interactive regression model") %&% tt_numeric_row(vdec.extract(rr$res_mini$vdec),perc=F,dec=1) +
    tt_text_row("mixture model") %&% tt_numeric_row(vdec.extract(rr$res_mixt$vdec),perc=F,dec=1) +
    tt_rule_bottom()

  tab = tt_tabularize(tt,header = "lrrrrr", pretty_rules=F)
  tt_save(tab,filename="inst/figures/tab-shimersmith-narrow_vdec.tex",stand_alone=F)

}



# =========================  EXTRA ANALYSIS =====================

test.liml.interactions <- function() {
  model_mini_bs         = archive.get("mini_bs",arch_static)
  M = laply(model_mini_bs,function(l) l$B1)
  M = M/M[,1]
  dd = data.table(melt(M[,2:10]))
  dd[, zval := ((value-1)/sd(value))^2,Var2]
  dd[, sum(zval) , Var1][,mean(V1)]

  qchisq(0.95,9)

  hist(dd$zval)
  dchisq
}

# =========================  DYNAMIC -- MINI =====================

# plot wages from the dynamic model
fig.dyn.unc <- function() {
  load("../figures/res/est-dynamic.dat")

  gp = wplot(res_unc$mode$A2s) + xlab("firm classes") + ylab("mean wage for each firm worker combination")
  ggsave("../figures/build/fig-dyn-mixt-unc-wage.pdf",gp,dpi=100,width = 8,height = 7)

  gp = pplot(res_unc$mode$pk0)+ xlab("firm classes") + ylab("worker type proportions")
  ggsave("../figures/build/fig-dyn-mixt-unc-qk.pdf",gp,dpi=100,width = 8,height = 7)
}

#' generates wages and sorting plots for the dynamic mini model
fig.dynamic.mini <- function() {
  model_mini         = archive.get("mini_model",arch_dyn)
  cstats             = archive.get("cstats",arch_dyn)

  model = model_mini; vs = c(-2,-1,0,1,2)
  rr = data.table(l=1:length(model$A1),Em = model$Em,Esd = c(model$Esd),
                  N = model$Ns,A1=model$A1,B1=model$B1,A2s=model$A2s,B2=model$B2)
  alpha_m  = rr[, wtd.mean(Em,N)]
  alpha_sd = sqrt(rr[, wtd.mean(Esd^2,N) + wtd.var(Em,N) ])

  rr2 = rr[, list( y1= (vs*alpha_sd + alpha_m)* B1+ A1 ,y2= (vs*alpha_sd + alpha_m)* B2+ A2s, k=1:length(vs)),l ]

  lmin = min(cstats[,m2-3*sd2])
  lmax = max(cstats[,m2+3*sd2])

  gp1 = ggplot(rr2,aes(x=l,y=y2)) + geom_line(aes(color=factor(k))) +
           theme_bw()   +theme(legend.position = "none") + scale_x_continuous("firm class k",breaks=1:10,labels=1:10,limits=c(1,10)) +scale_y_continuous("log-earnings") +
           coord_cartesian(ylim = c(lmin,lmax)) + geom_ribbon(aes(ymin=m2-2*sd2,ymax=m2+2*sd2,x=j1,y=m2),data=cstats,alpha=0.2)
  ggsave("../figures/build/fig-dynamic-mini-wages.pdf",plot = gp1,width = 6,height = 4)

  # this plots the proportions
  vs2 = c(-Inf,-1.5,-0.5,0.5,1.5,Inf)
  rr = data.table(expand.grid(k=1:10,s=1:length(vs)))
  rr[, pr := with(model,
                  pnorm(  ( vs2[s+1] - (Em[k]-alpha_m)/alpha_sd )/(Esd[k]/alpha_sd)  ) -
                    pnorm(  ( vs2[s]   - (Em[k]-alpha_m)/alpha_sd )/(Esd[k]/alpha_sd)  )) ,
     list(k,s)]
  gp2 = ggplot(rr,aes(x=k,y=pr,fill=factor(s))) + geom_bar(stat="identity") +
    theme_bw() +theme(legend.position = "none") + scale_x_continuous("firm class k",breaks=1:10,labels=1:10) +ylab("worker types")
  ggsave("../figures/build/fig-dynamic-mini-sorting.pdf",plot = gp2,width = 6,height = 4)

  pdf("../figures/build/fig-dynamic-mini.pdf",width = 8,height = 4)
  blm:::multiplot(gp1,gp2,cols = 2)
  dev.off()
}

fig.network <- function() {

  load("L:\\Tibo\\qtrdata\\tmp-2003-static.dat")
  jdata2 =jdata[,.N,list(f1,f2)]
  #fids = jdata2[,intersect(f1,f2)]
  #jdata2 = jdata2[f1 %in% fids][f2 %in% fids]

  jcount = jdata[,.N,list(f1,f2)]
  # make symetric
  jcount = rbind(jcount,jcount[,list(f1=f2,f2=f1,N=N)])
  jcount = jcount[,list(N=sum(N)),list(f1,f2)]

  # select larger firms
  fids = jcount[,sum(N),f1][V1>20,f1]

  jcount = jdata[f1 %in% fids][f2 %in% fids][,.N,list(f1,f2)]
  # make symetric
  jcount = rbind(jcount,jcount[,list(f1=f2,f2=f1,N=N)])
  jcount = jcount[,list(N=sum(N)),list(f1,f2)]

  jcount[,N:=N>=1]
  AA = acast(jcount,f1~f2,fill=0,drop=FALSE)
  net=graph.adjacency(AA,mode="undirected",weighted=NULL)

  plot(net,  vertex.size=3, vertex.label=NA,
       edge.arrow.size=0.3,
       rescale=TRUE)

  layout.by.attr <- function(graph, wc, cluster.strength=1,layout=layout.auto) {
    g <- graph.edgelist(get.edgelist(graph)) # create a lightweight copy of graph w/o the attributes.
    E(g)$weight <- 1

    attr <- cbind(id=1:vcount(g), val=wc)
    g <- g + vertices(unique(attr[,2])) + igraph::edges(unlist(t(attr)), weight=cluster.strength)

    l <- layout(g, weights=E(g)$weight)[1:vcount(graph),]
    return(l)
  }

  plot(net,  vertex.size=3, vertex.label=NA,
       edge.arrow.size=0.3,
       rescale=TRUE,layout=layout.by.attr(net, wc=1))

  ## voronoi plot
  load("L:\\Tibo\\qtrdata\\tmp-2003-static.dat")
  dd = sdata[,list(m1=mean(y1),m2=sd(y1),.N),f1][N>=10]
  dd[,plot(m1,m2)]

  ks = kmeansW(dd[,cbind(m1,m2)],centers = 10,weight = dd$N,iter.max = 1000,nstart = 500,scale())
  dd$centers = ks$cluster

  library(tripack)
  library(RColorBrewer)
  CL5 <- brewer.pal(10, "Set3")
  V <- voronoi.mosaic(ks$centers[,1],ks$centers[,2])
  P <- voronoi.polygons(V)
  plot(dd[N>50,cbind(m1,m2)],pch=19,xlim=c(9.6,11),ylim=c(0.15,0.5),
           xlab="wage mean",ylab="wage standard deviation",col=CL5[dd[N>50,centers]])
  #points(ks$centers[,1],ks$centers[,2],pch=3,cex=1.5,lwd=2)
  plot(V,add=TRUE)

  plot(dd[N>50,cbind(m1,m2)],pch=19,xlim=c(9.6,11),ylim=c(0.15,0.5),
       xlab="wage mean",ylab="wage standard deviation")


  dd = dd[]

  dd
  sdata[,x:=1][,y1_bu:=y1]
  sim = list(jdata=jdata,sdata=sdata)
  rm(sdata,jdata)

  # get clsuters
  grps  = res.load("m2-mixt-y2003-groups")
  sim   = grouping.append(sim,grps$best_cluster)
}

test.network.plot <- function() {
  m2stats = res.load("m2-stats-d2003")
  mstats = m2stats$mstats

  AA = acast(mstats,j1~j2,value=N)
  AA[AA<100]=0
  ig <- graph.adjacency(AA, mode="undirected", weighted=TRUE,diag = FALSE)
  plot(ig,edge.width=E(ig)$weight/100)

}
# =========================  STATIC -- PROBA =====================


plot.static.proba.gibbs <- function() {
  m2proba = res.load("m2-proba-gibbs-d2003-res")

  # we plot the variance of firm effect, and correlation
  m2proba$akm$vdec$model="akm"
  m2proba$akm$vdec$step=1:nrow(m2proba$akm$vdec)
  m2proba$blm$vdec$model="blm"
  m2proba$blm$vdec$step=1:nrow(m2proba$blm$vdec)
  rr = data.table(rbind(m2proba$akm$vdec,m2proba$blm$vdec))
  rr = rr[step<=165]

  rrm = data.table(melt(rr, id.vars = c("model","step","rsq") ))
  rrm[variable!="cor_kl",value := value * rsq]

  gp = ggplot(rrm,aes(x=step,y=value,color=factor(model),linetype=factor(model))) +
    geom_line() + theme_bw() + facet_wrap(~variable,scales = "free") +
    theme(legend.position = "none")
  ggsave("inst/figures/fig-proba-vdec.pdf",plot = gp,width = 6,height = 4)

  for (vv in c("var_l","var_k","cor_kl","cov_kl")) {
    gp = ggplot(rrm[variable==vv],aes(x=step,y=value,color=factor(model),linetype=factor(model))) +
      geom_line() + theme_bw() +
      theme(legend.position = "none") + ylab("") + xlab("")
    ggsave(sprintf("inst/figures/fig-proba-vdec-%s.pdf",vv),plot = gp,width = 4,height = 3)
  }


  # we plot the variance of firm effect, and correlation
  m2proba$akm$liks$model="akm"
  m2proba$akm$liks$step=1:nrow(m2proba$akm$liks)
  m2proba$blm$liks$model="blm"
  m2proba$blm$liks$step=1:nrow(m2proba$blm$liks)
  rr = data.table(rbind(m2proba$akm$liks,m2proba$blm$liks))
  rr = rr[outloop<=165]
  #rr[outloop==177 & model=="blm", lik := rr[outloop==178 & model=="blm"]$lik]
  #rr[outloop==177 & model=="blm", liks := rr[outloop==178 & model=="blm"]$liks]
  #rr[outloop==177 & model=="blm", likm := rr[outloop==178 & model=="blm"]$likm]

  rrm = melt(rr, id.vars = c("model","step"),measure.vars = c("lik","likm","liks","likg") )
  rrm[variable=="lik",variable:="total likelihood"][variable=="likm",variable:="movers likelihood"]
  rrm[variable=="liks",variable:="stayers likelihood"][variable=="likg",variable:="classification prior probabilities"]
  gp = ggplot(rrm,aes(x=step,y=value/1e6,color=factor(model),linetype=factor(model))) + geom_line() + theme_bw() + facet_wrap(~variable,scales = "free") + theme(legend.position = "none")
  ggsave("inst/figures/fig-proba-liks.pdf",plot = gp,width = 6,height = 4)

}

figure.hybrid <- function() {

  rr = data.table(res.load("m2-trace-d2003-withinRe"))
  rr[, akm_share := 0.32]

  gp = ggplot(rr,aes(x=nc,y= (cluster_var + omega_var_sub2)/y_var)) +
    geom_line() + theme_bw() + geom_line(aes(y=cluster_var/y_var),linetype=3) +
    geom_line(aes(y=akm_share),linetype=2) +
    ylim(0,0.4) + ylab("share of total variance") + xlab("number of firm classes") +
    annotate("text", x = 27, y = 0.335, label = "AKM firm variance") +
    annotate("text", x = 26, y = 0.01, label = "between-class variance") +
    annotate("text", x = 24, y = 0.07, label = "between-and-within-class variance")


  proj = list(path_figures = "inst/figures/")
  ggsave(file.path(proj$path_figures,"fig-within-re.pdf"),plot = gp,width = 6,height = 4)
}



test.reporting <- function() {
  rr = res.load("m2-blmhybrid")
  save(rr,file="../../smfe/inst/local-res/sweden-hybrid.dat")


}

tab.akm_discretization <- function() {
  rr = res.load("m2-akm-effect-of-discretization")

  tt =
    tt_text_row(c("K","L","$Var(\\psi)$","$Var(\\alpha)$","cov"))+
    tt_rule_mid() +
    tt_numeric_column(rr$nf,dec=0) %&%
    tt_numeric_column(rr$nk,dec=0) %&%
    tt_numeric_column(rr$vpsi,dec=3) %&%
    tt_numeric_column(rr$valpha,dec=3)%&%
    tt_numeric_column(rr$cov,dec=3)

  tab = tt_tabularize(tt,header = "llrrr", pretty_rules=T)
  tt_save(tab,filename="inst/figures/tab-akm-discretization-effect.tex",stand_alone=F)
}

tab.dirty_iteration <- function() {
  rr = data.table(res.load("m2-mixt-d2003-main-reclassify-allresults"))

  lseq = paste(c(0,1,2,3,5,10,15,20))

  tt_user = function(key,lseq,title) {
    with(rr[start==key][.id %in% lseq],{
      tt_text_row(title,cspan = 6,center = "c") +
      tt_spacer_row(-2) +
      tt_rule_mid() +
      tt_text_column(.id) %&%
        tt_numeric_column(100*var_k*rsq,dec=2,percent=F) %&%
        tt_numeric_column(100*var_l*rsq,dec=2,percent=F) %&%
        tt_numeric_column(100*cov_kl*rsq,dec=2,percent=F) %&%
        tt_numeric_column(100*(1-rsq),dec=2,percent=F) %&%
        tt_numeric_column(100*cor_kl,dec=2,percent=F)
    })
  }

  tt =
    tt_text_row(c("iter.","$\\frac{Var(\\alpha)}{Var(y)}$",
                  "$\\frac{Var(\\psi)}{Var(y)}$",
                  "$\\frac{2 Cov(\\alpha,\\psi)}{Var(y)}$",
                  "$\\frac{Var(\\varepsilon)}{Var(y)}$",
                  "$Corr(\\alpha,\\psi)$")) +
    tt_rule_mid() +
    tt_spacer_row(5) +
    #tt_user("mean",lseq,"Starting using means") +
    #tt_spacer_row(5) +
    tt_user("va",lseq,"Starting using value added") +
    tt_spacer_row(5) +
    tt_user("prank",lseq,"Starting using poaching rank") +
    tt_spacer_row(5) +
    tt_user("akm",lseq,"Starting using AKM") +
    tt_spacer_row(5) +
    tt_user("mratio",lseq,"Starting using share of movers") +
    tt_spacer_row(5) +
    tt_user("resid",lseq,"Using residual earnings")

  tab = tt_tabularize(tt,header = "rrrrrr", pretty_rules=T)
  tt_save(tab,filename="inst/figures/tab-dirty-iteration.tex",stand_alone=F)

  tt =
    tt_rule_top() +
    tt_text_row(c("iter.","$\\frac{Var(\\alpha)}{Var(y)}$",
                  "$\\frac{Var(\\psi)}{Var(y)}$",
                  "$\\frac{2 Cov(\\alpha,\\psi)}{Var(y)}$",
                  "$\\frac{Var(\\varepsilon)}{Var(y)}$",
                  "$Corr(\\alpha,\\psi)$")) +
    tt_rule_mid() +
    #tt_user("mean",lseq,"Starting using means") +
    #tt_spacer_row(5) +
    tt_user("dynanmic",lseq,"Dynamic model") +
    tt_rule_bottom()

  tab = tt_tabularize(tt,header = "rrrrrr", pretty_rules=F)
  tt_save(tab,filename="inst/figures/tab-dirty-iteration-dyn.tex",stand_alone=F)

}

