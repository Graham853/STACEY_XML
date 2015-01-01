


add.distribution <- function(A) {
  # TODO id as function
  add.opennode("distribution", attrs=c(id="posterior", spec="util.CompoundDistribution")) 
  add.coalescent.distribution(A)
  add.prior.distribution(A)
  add.likelihood.distribution(A)
  add.closetag()
}



add.coalescent.distribution <- function(A) {
  add.bigcomment("the STACEY pops-integrated-out coalescent")
  attrs <- c(id=smcCoalescentID(), spec="stacey.PIOMSCoalescentDistribution", 
             tree=IDtoREF(smcTreeID()), taxonset=IDtoREF(taxonSetOfSetsID()))
  add.opennode("distribution", attrs=attrs)
  for (g in 1:nof.alignments()) {
    # TODO ploidy
    attrs <- c(id=geneTreeCoalFactorID(g), spec="stacey.GtreeAndCoalFactor", 
               tree=IDtoREF(geneTreeID(g)), Ploidy="2.0")
    add.node("geneTree", attrs=attrs)
  }
  
  smccoal <- get.smc.coalescent()
  w <- smccoal$invgammamix$weights
  a <- smccoal$invgammamix$alphas
  b <- smccoal$invgammamix$betas
  for (i in 1:length(w)) {
    attrs <- c(id=InvGammaComponentID(i), spec="stacey.InverseGammaComponent", 
               weight=w[i], alpha=a[i], beta=b[i])
    add.node("popPriorInvGamma", attrs=attrs)
  }
  add.node(popSFID(), idref=IDtoREF(popSFID())) 
  add.closetag()  
}


add.likelihood.distribution <- function(A) {
  add.bigcomment("Gene tree likelihoods")
  add.opennode("distribution", attrs=c(id="likelihood", spec="util.CompoundDistribution"))
  for (g in 1:nof.alignments()) {
    add.gtree.lhood(g)
  } 
  add.closetag()
}


add.gtree.lhood <- function(g) {
  add.comment(paste0("tree likelihood for gene ", g))
  attrs <- c(id=geneTreeLhoodID(g), data=IDtoREF(alignmentID(g)), 
             spec="TreeLikelihood", tree="TreeLikelihood")
  add.opennode("distribution", attrs=attrs)
  add.comment(paste0("site model for gene ", g))
  add.sitemodel(g)
  add.comment(paste0("branch model for gene ", g))
  add.clockmodel(g)
  add.closetag()
}


  
add.sitemodel <- function(g) { 
  #TODO this doesn't refer properly to site model when two genes share a site model
  
  sitehets <- get.sitehets()
  substs <- get.substs()
  siteID <- paste0("SiteModel.", g)
  muID <- paste0("mutationRate.", g)
  gamID <- paste0("siteHetGamma.", g)
  invID <- paste0("proportionInvariant.", g)
  hkyID <- paste0("hky.", g)
  freqsID <- paste0("estimatedFreqs.", g)
  
  add.opennode("siteModel", attrs=c(id=siteID, spec="SiteModel"))
  add.node("parameter", attrs=c(id=muID, name="mutationRate", estimate="false"), .children="1.0")
  add.node("parameter", attrs=c(id=gamID, name="shape", estimate="false"), .children="1.0")
  attrs <- c(id=invID, name="proportionInvariant", estimate="false", lower="0.0", upper="1.0")
  add.node("parameter", attrs=attrs, .children="0.0")
  add.closetag()
  
  add.opennode("substModel", attrs=c(id=hkyID, kappa=IDtoREF(kappaID(g)), spec="HKY"))
  attrs <- c(id=freqsID, frequencies=IDtoREF(frequenciesID(g)), spec="Frequencies")
  add.node("frequencies", attrs=attrs)
  add.closetag()
}




add.clockmodel <- function(g) {
  #TODO this doesn't refer properly to clock model when two genes share a clock
  scID <- paste0("StrictClock.", g)
  if (g == 1) {
    # in g==1 case make a new fixed parameter
    attrs <- c(id=scID, spec="beast.evolution.branchratemodel.StrictClockModel")
    add.opennode("branchRateModel", attrs=attrs)
    attrs <- c(id=clockRateID(g), estimate="false", name="clock.rate")
    add.node("parameter", attrs=attrs, .children="1.0")
    add.closetag()
  } else {
    # else refer to existing clock rate
    attrs <- c(id=scID, clock.rate=IDtoREF(clockRateID(g)), spec="beast.evolution.branchratemodel.StrictClockModel")
    add.opennode("branchRateModel", attrs=attrs)
    add.closetag()
  }
}




add.prior.distribution <- function(A) {
  
  add.bigcomment("prior")
  add.opennode("distribution", attrs=c(id="prior", spec="util.CompoundDistribution"))
  add.comment("BirthDeathCollapse prior for smcTree") 
  
  attrs <- c(id=bdcModelID(), spec="stacey.BirthDeathCollapseModel",
             collapseHeight=get.bdc.model()$eps)
  add.opennode("distribution", attrs=attrs)
  add.node("tree", attrs=c(id=smcTreeID()))
  add.node("parameter", attrs=c(id=bdcGrowthID(), name="birthDiffRate"))
  add.node("parameter", attrs=c(id=bdcRelDeathID(), name="relativeDeathRate"))
  add.node("parameter", attrs=c(id=bdcCollapseWtID(), name="collapseWeight"))
  add.node("parameter", attrs=c(id=bdcOriginHtID(), name="collapseWeight"))    
  add.closetag()
  
  
  
  add.comment("Priors for BirthDeathCollapse hyperparameters") 
  add.1Dprior(get.bdc.model()$growthrate, IDtoREF(bdcGrowthID()))   
  add.1Dprior(get.bdc.model()$reldeath,   IDtoREF(bdcRelDeathID()))
  add.1Dprior(get.bdc.model()$w,          IDtoREF(bdcCollapseWtID()))
   
  TheSTACEYxmlTree$addComment("Hyper-prior for population scale factor") 
  add.1Dprior(get.smc.coalescent()$popSF, IDtoREF(popSFID()))

  add.comment("Priors for relative clock rates of gene trees (after first one)")
  clocks <- get.clocks()
  if (length(clocks) > 1) {
    for (g in 2:length(clocks)) {
      add.1Dprior(clocks[[g]]$clock$rate, IDtoREF(clockRateID(g)))   
    }    
  }
  
  add.comment("Priors for substitution models")
  substs <- get.substs()
  for (g in 1:length(substs)) {
    if (substs[[g]]$model$kind == "HKY") {
      add.1Dprior(substs[[g]]$model$kappa, IDtoREF(kappaID(g)))
      if (substs[[g]]$model$freqs$kind == "Dirichlet"  &&
            dirichlet.is.uniform(substs[[g]]$model$freqs$mean, substs[[g]]$model$freqs$alpha)) {
        # nothing to do
      }
    } else {
      stop("Only prior for frequencies implemented is uniform Dirichlet")
    }
  }  
  
add.closetag()
}


dirichlet.is.uniform <- function(mean, alpha) {
  if (alpha == 1) {
    if (length(mean) > 0) {
      ok = TRUE
      for (i in 1:length(mean)) {
        ok <- ok && isTRUE(all.equal(mean[1], mean[i], tolerance=1e-12))
      }
      ok
    } else {
      FALSE
    }
  } else {
    FALSE
  }
}














