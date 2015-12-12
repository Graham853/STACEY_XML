


add.distribution <- function() {
  open.xmlnode("distribution", attrs=c(id=posteriorID(), spec="util.CompoundDistribution")) 
  add.coalescent.distribution()
  add.prior.distribution()
  add.likelihood.distribution()
  close.xmlnode()
}



add.coalescent.distribution <- function() {
  add.bigcomment("the STACEY pops-integrated-out coalescent")
  attrs <- c(id=smcCoalescentID(), spec="stacey.PIOMSCoalescentDistribution", 
             tree=IDtoREF(smcTreeID()), taxonset=IDtoREF(taxonSetOfSetsID()))
  open.xmlnode("distribution", attrs=attrs)
  gtrees <- get.GTrees()
  for (u in 1:length(gtrees)) {
    attrs <- c(id=geneTreeCoalFactorID.u(u), spec="stacey.GtreeAndCoalFactor", 
               tree=IDtoREF(geneTreeID.u(u)), Ploidy=gtrees[[u]]$ploidy)
    add.xmlnode("geneTree", attrs=attrs)
  }
  
  smccoal <- get.smc.coalescent()
  popSF <- get.child(smccoal, "popSF") 
  if (is.estimated(popSF)) {
    add.xmlnode("popPriorScale", attrs=c(idref=popSFID()))  
  } else {
    add.fixed.parameter(popSFID(), "popPriorScale", popSF$value)
  }
  ivgm <- get.child(smccoal, "invgammamix") 
  w <- get.childvalue(ivgm, "weights")
  a <- get.childvalue(ivgm, "alphas")
  b <- get.childvalue(ivgm, "betas")
  for (i in 1:length(w)) {
    attrs <- c(id=InvGammaComponentID(i), spec="stacey.InverseGammaComponent", 
               weight=w[i], alpha=a[i], beta=b[i])
    add.xmlnode("popPriorInvGamma", attrs=attrs)
  }
  
  close.xmlnode()  
}





add.likelihood.distribution <- function() {
  add.bigcomment("Tree likelihoods for loci")
  open.xmlnode("distribution", attrs=c(id=likelihoodID(), spec="util.CompoundDistribution"))
  for (a in 1:get.nof.alignments()) {
    add.partition.lhood(a)
  } 
  close.xmlnode()
}


add.partition.lhood <- function(a) {
  add.comment(paste0("tree likelihood for partition ", a))
  # When site models are linked, the first occurence is added as an element,
  # and others refer to it   
  siteMs <- get.gtree.SiteMs.with.first.partition()
  siteM.isfirst <- FALSE
  for (u in 1:length(siteMs)) {
    if (siteMs[[u]]$first.partition == a) {
      siteM.isfirst <- TRUE
    }
  }
  
  if (siteM.isfirst) {
    attrs <- c(id=geneTreeLhoodID(a), data=IDtoREF(partitiondataID.a(a)), 
               spec="TreeLikelihood", tree=IDtoREF(geneTreeID.a(a)))
  } else {
    attrs <- c(id=geneTreeLhoodID(a), data=IDtoREF(partitiondataID.a(a)), 
               spec="TreeLikelihood", tree=IDtoREF(geneTreeID.a(a)),
               siteModel=IDtoREF(sitemodelID.a(a)))
  }
  open.xmlnode("distribution", attrs=attrs)
  add.comment(paste0("site model for partition ", a))
  if (siteM.isfirst) {
    add.sitemodel(a)
  }
  add.comment(paste0("branch model for partition ", a))
  add.branchratemodel(a)
  close.xmlnode()
}


  
add.sitemodel <- function(a) {  
	substM <- get.gtree.SubstM(a)
	sitehet <- get.gtree.SiteHet(a)
  open.xmlnode("siteModel", attrs=c(id=sitemodelID.a(a), spec="SiteModel"))
  if (sitehet$value == "None") {
    add.xmlnode.children("parameter", attrs=c(id=sitehetMuID.a(a), name="mutationRate", estimate="false"), children="1.0")
    add.xmlnode.children("parameter", attrs=c(id=sitehetGammaShapeID.a(a), name="shape", estimate="false"), children="1.0")
    attrs <- c(id=sitehetPropInvID.a(a), name="proportionInvariant", estimate="false", lower="0.0", upper="1.0")
    add.xmlnode.children("parameter", attrs=attrs, children="0.0")
  } else {
    # TODO other gamma counts, other site het models (if any)
  }
	model <- get.child(substM, "model")
  if (model$kind == "HKY") {
  	kappa <- get.childvalue(model, "kappa")
  	freqs <- get.childvalue(model, "freqs")
    if (is.numeric(kappa)) {
      if (is.character(freqs)) {
        # fixed kappa and freqs
        open.xmlnode("substModel", attrs=c(id=substsmodelID.a(a), spec="HKY"))
        kappa.attrs <- c(name="kappa", lower="0.0", estimate="false")
        add.xmlnode.children("parameter", attrs=kappa.attrs, children=kappa)
        # TODO need freqs$kind here(?)
        freqs.attrs <- c(id=empiricalFreqsID.a(a), spec="Frequencies", data=IDtoREF(partitiondataID.a(a)))
        add.xmlnode("frequencies", attrs=freqs.attrs)  
        close.xmlnode()
      } else {
        # TODO fixed kappa, estimated freqs
        stop("Unimplemented option in HKY: fixed kappa, estimated freqs")
      }
    } else {
      if (is.character(freqs)) {
        # TODO estimated kappa, fixed freqs
        stop("Unimplemented option in HKY: estimated kappa, fixed freqs")
      } else {
        # estimated kappa and freqs
        open.xmlnode("substModel", attrs=c(id=substsmodelID.a(a), kappa=IDtoREF(kappaID.a(a)), spec="HKY"))
        attrs <- c(id=frequenciesModelID.a(a), frequencies=IDtoREF(frequenciesParamID.a(a)), spec="Frequencies")
        add.xmlnode("frequencies", attrs=attrs)
        close.xmlnode()
      }
    }
         
  } else {
    # TODO other subst models
    stop("Unimplemented substitution model ", siteM$subst$model$kind)
  }
  
  close.xmlnode() 
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
    open.xmlnode("branchRateModel", attrs=attrs)
    attrs <- c(id=clockRateID.a(a), estimate="false", name="clock.rate")
    add.xmlnode.children("parameter", attrs=attrs, children="1.0")
    close.xmlnode()
  } else {
    # else refer to existing clock rate
    attrs <- c(id=scID, clock.rate=IDtoREF(clockRateID.a(a)), spec="beast.evolution.branchratemodel.StrictClockModel")
    open.xmlnode("branchRateModel", attrs=attrs)
    close.xmlnode()
  }
}




add.prior.distribution <- function() {
  
  add.bigcomment("prior")
  open.xmlnode("distribution", attrs=c(id="prior", spec="util.CompoundDistribution"))
  add.comment("BirthDeathCollapse prior for smcTree") 
  
  attrs <- c(id=bdcModelID(), spec="stacey.BirthDeathCollapseModel",
             collapseHeight=get.childvalue(get.bdc.model(), "eps"))
  open.xmlnode("distribution", attrs=attrs)
  add.xmlnode("tree", attrs=c(idref=smcTreeID()))
  
  add.comment("Parameters (growth=lambda-mu, reldeath=mu/lambda, collapse weight=w)") 
  add.comment("for the BirthDeathCollapse model") 
  bdc <- get.bdc.model()
  
  bdc.g <- get.child(bdc, "growthrate")
  if (is.estimated(bdc.g)) {
    add.xmlnode("parameter", attrs=c(idref=bdcGrowthID(), name="birthDiffRate")) 
  } else {
    add.fixed.parameter(bdcGrowthID(), "birthDiffRate", bdc.g$value)
  }
  bdc.rd <- get.child(bdc, "reldeath")
  if (is.estimated(bdc.rd)) {
    add.xmlnode("parameter", attrs=c(idref=bdcRelDeathID(), name="relativeDeathRate"))
  } else {
    add.fixed.parameter(bdcRelDeathID(), "relativeDeathRate", bdc.rd$value)
  }
  bdc.w <- get.child(bdc, "w")
  if (is.estimated(bdc.w)) {
    add.xmlnode("parameter", attrs=c(idref=bdcCollapseWtID(), name="collapseWeight")) 
  } else {
    add.fixed.parameter(bdcCollapseWtID(), "collapseWeight", bdc.w$value) 
  }
  add.comment("Origin height. Always estimated using a built-in prior.")
  add.xmlnode("parameter", attrs=c(idref=bdcOriginHtID(), name="originHeight"))    
  close.xmlnode()
  add.comment("Priors for the BirthDeathCollapse model parameters (if they are estimated)") 
  if (is.estimated(bdc.g)) {
    add.1Dprior(bdc.g,   IDtoREF(bdcGrowthID()))
  }
  if (is.estimated(bdc.rd)) {
    add.1Dprior(bdc.rd,   IDtoREF(bdcRelDeathID()))
  }
  if (is.estimated(bdc.w)) {
    add.1Dprior(bdc.w,          IDtoREF(bdcCollapseWtID()))
  }
  
   
  popSF <- get.child(get.smc.coalescent(), "popSF")
  if (is.estimated(popSF)) {
    add.comment("Prior for population scale factor") 
    add.1Dprior(popSF, IDtoREF(popSFID()))
  } 

  add.comment("Priors for relative clock rates of partitions (if estimated)")
  clocks <- get.Clocks()
  for (u in 1:length(clocks)) {
  	prm <- get.child(clocks[[u]], "PartitionRateM")
  	clkrate <- get.child(prm, "rate")
    if (is.estimated(clkrate)) {
      add.1Dprior(clkrate, IDtoREF(clockRateID.u(u)))
    } 
  }
  
  add.comment("Priors for substitution models")
  siteMs <- get.SiteMs()
  for (u in 1:length(siteMs)) {
  	substM <- get.child(siteMs[[u]], "SubstM")
  	model <- get.child(substM, "model")
    if (model$kind == "HKY") {
    	kappa <- get.child(model, "kappa")
    	freqs <- get.child(model, "freqs")
      if (is.numeric(kappa$value)) {
        #  nothing to do
      } else {
        add.1Dprior(kappa, IDtoREF(kappaID.u(u)))
      }
      if (is.character(freqs$value)) {
        #  nothing to do
      } else {
        if (freqs$kind == "UniformUnitSimplex") {
          # still nothing to do
        } else {
          stop("Only prior for frequencies implemented is UniformUnitSimplex")
        }
      }
    } 
  }  
  
close.xmlnode()
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














