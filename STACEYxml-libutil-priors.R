


add.1Dprior <- function(prior, xREF) {
  open.xmlnode("prior", attrs=c(name="distribution", x=xREF))
  add.ParametricDistribution(prior)
  close.xmlnode()
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



add.fixed.parameter <- function(id, name, value) {
  add.xmlnode.children("parameter", attrs=c(id=id, name=name, estimate="false"), children=value) 
}


# <map name="Beta">beast.math.distributions.Beta</map>
#   <map name="Exponential">beast.math.distributions.Exponential</map>
#   <map name="InverseGamma">beast.math.distributions.InverseGamma</map>
#   <map name="LogNormal">beast.math.distributions.LogNormalDistributionModel</map>
#   <map name="Gamma">beast.math.distributions.Gamma</map>
#   <map name="Uniform">beast.math.distributions.Uniform</map>
#   <map name="LaplaceDistribution">beast.math.distributions.LaplaceDistribution</map>
#   <map name="OneOnX">beast.math.distributions.OneOnX</map>
#   <map name="Normal">beast.math.distributions.Normal</map>

add.lnorm.prior <- function(prior) {
  open.xmlnode("LogNormal", attrs=c(id=paste0("LogNormal.", prior$id), name="distr"))
  add.fixed.parameter(paste0("meanlog.", prior$id), "M", get.childvalue(prior, "meanlog"))
  add.fixed.parameter(paste0("sdlog.", prior$id), "S", get.childvalue(prior, "sdlog"))
  close.xmlnode()
}


add.beta.prior <- function(prior) {
  open.xmlnode("Beta", attrs=c(id=paste0("Beta.", prior$id), name="distr"))
  add.fixed.parameter(paste0("alpha.", prior$id), "alpha", get.childvalue(prior, "a"))
  add.fixed.parameter(paste0("beta.", prior$id), "beta", get.childvalue(prior, "b"))
  close.xmlnode()
}



add.exp.prior <- function(prior) {
  open.xmlnode("Exponential", attrs=c(id=paste0("Exponential.", prior$id), name="distr"))
  #TODO mean or scale
  close.xmlnode()
}


# TODO something funny about uniform limits
add.unif.prior <- function(prior) {
  open.xmlnode("Uniform", attrs=c(id=paste0("Uniform.", prior$id), name="distr"))
  add.fixed.parameter(paste0("lower.", prior$id), "lower", get.childvalue(prior, "lower"))
  add.fixed.parameter(paste0("upper.", prior$id), "upper", get.childvalue(prior, "upper"))
  close.xmlnode()
}


