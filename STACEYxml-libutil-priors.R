


add.1Dprior <- function(prior, xREF) {
  add.opennode("prior", attrs=c(name="distribution", x=xREF))
  add.ParametricDistribution(prior)
  add.closetag()
}


add.ParametricDistribution <- function(prior) {
  if (identical(prior$kind, "LogNorm")) {
    add.lnorm.prior(prior)
  } else if (identical(prior$kind, "Beta")) {
    add.beta.prior(prior)
  } else if (identical(prior$kind, "Exponential")) {
    stop("not implemented in make.ParametricDistribution()")
    add.exp.prior(prior)
  } else if (identical(prior$kind, "Uniform")) {
    add.unif.prior(prior)
  } else {
    stop(paste0("unknown prior (kind=", prior$kind, ", id=", prior$id, ") in make.ParametricDistribution()"))
  }  
}




add.lnorm.prior <- function(prior) {
  add.opennode("LogNormal", attrs=c(id=paste0("LogNormal.", prior$id), name="distr"))
  attrs <- c(id=paste0("meanlog.", prior$id), name="M", estimate="false")
  add.node("parameter", attrs=attrs, .children=prior$meanlog)
  attrs <- c(id=paste0("sdlog.", prior$id),   name="S", estimate="false")
  add.node("parameter", attrs=attrs, .children=prior$sdlog)
  add.closetag()
}


add.beta.prior <- function(prior) {
  add.opennode("Beta", attrs=c(id=paste0("Beta.", prior$id), name="distr"))
  attrs <- c(id=paste0("alpha.", prior$id), name="alpha", estimate="false")
  add.nodee("parameter", attrs=attrs, .children=prior$alpha)
  attrs <- c(id=paste0("beta.", prior$id), name="beta", estimate="false")
  add.node("parameter", attrs=attrs, .children=prior$beta)
  add.closetag()
}



add.exp.prior <- function(prior) {

  add.opennode("Exponential", attrs=c(id=paste0("Exponential.", prior$id), name="distr"))
  #TODO mean or scale
  add.closetag()
  
}



add.unif.prior <- function(prior) {
  add.opennode("Uniform", attrs=c(id=paste0("Uniform.", prior$id), name="distr"))
  attrs <- c(id=paste0("lower.", prior$id), name="lower", estimate="false")
  add.node("parameter", attrs=attrs, .children=prior$lower)
  attrs <- c(id=paste0("upper.", prior$id), name="upper", estimate="false")
  add.node("parameter", attrs=attrs, .children=prior$upper)
  add.closetag()
}

