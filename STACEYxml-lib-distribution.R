


add.distribution <- function(A) {
  # TODO id as function
  add.opennode("distribution", attrs=c(id=posteriorID(), spec="util.CompoundDistribution")) 
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
  gtrees <- get.gtrees()
  for (g in 1:length(gtrees)) {
    # TODO ploidy
    attrs <- c(id=geneTreeCoalFactorID.g(g), spec="stacey.GtreeAndCoalFactor", 
               tree=IDtoREF(geneTreeID.g(g)), Ploidy="2.0")
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
  add.node("popPriorScale", attrs=c(idref=popSFID())) 
  add.closetag()  
}


add.likelihood.distribution <- function(A) {
  add.bigcomment("Tree likelihoods for loci")
  add.opennode("distribution", attrs=c(id=likelihoodID(), spec="util.CompoundDistribution"))
  for (a in 1:nof.alignments()) {
    add.partition.lhood(A, a)
  } 
  add.closetag()
}


add.partition.lhood <- function(A, a) {

  
  add.comment(paste0("tree likelihood for partition ", a))
  attrs <- c(id=geneTreeLhoodID(a), data=IDtoREF(geneTreeID.a(a)), 
             spec="TreeLikelihood", tree="TreeLikelihood")
  add.opennode("distribution", attrs=attrs)
  add.comment(paste0("site model for partition ", a))
  add.sitemodel(A, a)
  add.comment(paste0("branch model for partition ", a))
  add.branchratemodel(a)
  add.closetag()
}


  
add.sitemodel <- function(A, a) {  
  # When site models are linked, the first occurence is added as an element,
  # and others refer to it 
  
  # <distribution data="@U2SQg5t10r1-2" id="treeLikelihood.U2SQg5t10r1-2" siteModel="@SiteModel.s:U2SQg5t10r1-1" spec="TreeLikelihood" tree="@Tree.t:U2SQg5t10r1-2">
  #  <branchRateModel clock.rate="@clockRate.c:U2SQg5t10r1-2" id="StrictClock.c:U2SQg5t10r1-2" spec="beast.evolution.branchratemodel.StrictClockModel"/>
  # </distribution>
  
  siteMs <- get.siteMs()
  first <- FALSE
  for (u in 1:length(siteMs)) {
    if (siteMs[[u]]$first.partition == a) {
      first <- TRUE
    }
  }
  siteM <- A$alignment.table$alignments[[a]]$siteM
  stopifnot(is.character(siteM$sitehet$model))
  if (first > 0) {
    add.opennode("siteModel", attrs=c(id=sitemodelID.a(a), spec="SiteModel"))
    if (siteM$sitehet$model == "None") {
      add.node("parameter", attrs=c(id=sitehetMuID.a(a), name="mutationRate", estimate="false"), .children="1.0")
      add.node("parameter", attrs=c(id=sitehetGammaShapeID.a(a), name="shape", estimate="false"), .children="1.0")
      attrs <- c(id=sitehetPropInvID.a(a), name="proportionInvariant", estimate="false", lower="0.0", upper="1.0")
      add.node("parameter", attrs=attrs, .children="0.0")
      
      add.opennode("substModel", attrs=c(id=substsmodelID.a(a), kappa=IDtoREF(kappaID.a(a)), spec="HKY"))
      attrs <- c(id=frequenciesModelID.a(a), frequencies=IDtoREF(frequenciesParamID.a(a)), spec="Frequencies")
      add.node("frequencies", attrs=attrs)
      add.closetag()
      add.closetag()      
    } else {
      # TODO
    }

  } else {
    
  }
}

# From Beauti
# here, trees 3,4,5 are linked.
# clock rates (partition rates) are all unlinked.
# site models are omitted from XML.
# clock model (branch rate model) is strict except for linked tree, which is ucln.
# The ucln not 'normalized' (normalized is just a flag to the ucln)
# clock rate 2,4,5 are estimated. (they are relative to 1 so no 1. 3 is estimated as part of ucln)

# <distribution data="@U2SQg5t10r1-1" id="treeLikelihood.U2SQg5t10r1-1" spec="TreeLikelihood" tree="@Tree.t:U2SQg5t10r1-1">
# <branchRateModel id="StrictClock.c:U2SQg5t10r1-1" spec="beast.evolution.branchratemodel.StrictClockModel">
#   <parameter estimate="false" id="clockRate.c:U2SQg5t10r1-1" name="clock.rate">1.0</parameter>
# </branchRateModel>
# </distribution>
# 
# <distribution data="@U2SQg5t10r1-2" id="treeLikelihood.U2SQg5t10r1-2" spec="TreeLikelihood" tree="@Tree.t:U2SQg5t10r1-2">
# <branchRateModel clock.rate="@clockRate.c:U2SQg5t10r1-2" id="StrictClock.c:U2SQg5t10r1-2" spec="beast.evolution.branchratemodel.StrictClockModel"/>
# </distribution>
# 
# <distribution data="@U2SQg5t10r1-3" id="treeLikelihood.U2SQg5t10r1-3" spec="TreeLikelihood" tree="@Tree.t:U2SQg5t10r1-3">
# <branchRateModel clock.rate="@ucldMean.c:U2SQg5t10r1-3" id="RelaxedClock.c:U2SQg5t10r1-3" rateCategories="@rateCategories.c:U2SQg5t10r1-3" spec="beast.evolution.branchratemodel.UCRelaxedClockModel" tree="@Tree.t:U2SQg5t10r1-3">
#   <LogNormal S="@ucldStdev.c:U2SQg5t10r1-3" id="LogNormalDistributionModel.c:U2SQg5t10r1-3" meanInRealSpace="true" name="distr">
#     <parameter estimate="false" id="RealParameter.01" lower="0.0" name="M" upper="1.0">1.0</parameter>
#   </LogNormal>
# </branchRateModel>
# </distribution>
# 
# <distribution id="treeLikelihood.U2SQg5t10r1-4" spec="TreeLikelihood" tree="@Tree.t:U2SQg5t10r1-3">
# <data idref="U2SQg5t10r1-4"/>
# <branchRateModel clock.rate="@clockRate.c:U2SQg5t10r1-4" id="StrictClock.c:U2SQg5t10r1-4" spec="beast.evolution.branchratemodel.StrictClockModel"/>
# </distribution>
# 
# <distribution id="treeLikelihood.U2SQg5t10r1-5" spec="TreeLikelihood" tree="@Tree.t:U2SQg5t10r1-3">
# <data idref="U2SQg5t10r1-5"/>
# <siteModel id="SiteModel.s:U2SQg5t10r1-5" spec="SiteModel">
# <branchRateModel clock.rate="@clockRate.c:U2SQg5t10r1-5" id="StrictClock.c:U2SQg5t10r1-5" spec="beast.evolution.branchratemodel.StrictClockModel"/>
# </distribution>



add.branchratemodel <- function(a) {
  # TODO this only does strict clock
  # When trees are linked, all but the first get a strict clock model.
  # So always get a branch rate model.
  # clock rates (partition rates) 
  scID <- paste0("StrictClock.", a)
  if (a == 1) {
    # in a==1 case make a new fixed parameter
    attrs <- c(id=scID, spec="beast.evolution.branchratemodel.StrictClockModel")
    add.opennode("branchRateModel", attrs=attrs)
    attrs <- c(id=clockRateID.a(a), estimate="false", name="clock.rate")
    add.node("parameter", attrs=attrs, .children="1.0")
    add.closetag()
  } else {
    # else refer to existing clock rate
    attrs <- c(id=scID, clock.rate=IDtoREF(clockRateID.a(a)), spec="beast.evolution.branchratemodel.StrictClockModel")
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

  add.comment("Priors for relative clock rates of partition trees (after first one)")
  clocks <- get.clocks()
  if (length(clocks) > 1) {
    for (c in 2:length(clocks)) {
      add.1Dprior(clocks[[c]]$clock$rate, IDtoREF(clockRateID.c(c)))   
    }    
  }
  
  add.comment("Priors for substitution models")
  siteMs <- get.siteMs()
  for (u in 1:length(siteMs)) {
    subst <- siteMs[[u]]$subst$model
    if (subst$kind == "HKY") {
      add.1Dprior(subst$kappa, IDtoREF(kappaID.u(u)))
      if (subst$freqs$kind == "UniformUnitSimplex") {
        # nothing to do
      }
    } else {
      stop("Only prior for frequencies implemented is UniformUnitSimplex")
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














