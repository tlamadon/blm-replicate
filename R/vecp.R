# serialization of parameters
# set of functions that allows to set and extract 
# parameters from a vector.

vecp <- function() {
  pp = list(end=0,p=list())
  class(pp) = "vecp"
  return(pp)
}

print.vecp <- function(pp) {
  for (p in pp$p) {
    cat(sprintf("%s %s ( %i values )\n",p$name,paste(p$dim,collapse="x"),p$length))
  }
}

v_param <- function(name,dim=1,ni=NULL,nv=NULL) {
  rr        = list(name=name,dim=dim,ni=ni,nv=nv)
  rr$length = cumprod(dim)[length(dim)] - length(ni)
  class(rr) <- "vecparam"
  return(rr)
}

`+.vecp` <- function(pp, rr){
  if (inherits(rr,"vecparam")) {
    # append to pp
    rr$start  = pp$end
    pp$end    = pp$end + rr$length
    pp$p[[ rr$name ]] = rr
  }
  return(pp)
}

setp <- function(pp,v,name) {
  
}

getp <- function(pp,v,name) {
  p = pp$p[[name]]
  I = (p$start+1):(p$start+p$length)
  R = v[I]
  dim(R) <- p$dim
  return(R)
}

getis <- function(pp,name) {
  p = pp$p[[name]]
  return(p$start:(p$start+p$length-1))
}


geti <- function(pp,name,coord) {
  p = pp$p[[name]]
  ll = p$start + coord[length(p$dim)] -1
  if (length(p$dim)>1) for (i in (length(p$dim)-1):1) {
    ll = coord[i] -1 + p$dim[i]*ll
  }
  return(ll+1)
}

test <- function() {
  
  vp <- vecp() + 
         v_param("A",c(5,5),c(1,2),c(1.0,2.0)) +
         v_param("B",c(6),c(2),c(3.0))
   
}