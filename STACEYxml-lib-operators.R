
add.operators <- function(A) {
  nofmcs <- length(get.minclusters())
  add.comment("Operators for both smcTree and locus trees.")
  add.nodeReheight(18)
  if (nofmcs > 2) {
    add.nodesNudge(6)
    add.coordinatedPruneRegraft(6)
  }
  if (nofmcs > 3) {
    add.focusedScaler(6)
  }
  
  add.comment("stretch and squeeze everything") 
  add.upGrowthClocks.downPopsHeights(2) 
  
  add.comment("Scalers for the birth-death-collapse model (hyper)parameters") 
  if (!get.bdcm.growthrate.fixed()) 
    { add.general.scaleOperator(bdcGrowthScalerID(), bdcGrowthID(), 0.75, 1) }
  if (!get.bdcm.reldeath.fixed())
    { add.general.scaleOperator(bdcReldeathScalerID(), bdcRelDeathID(), 0.75, 1) }
  if (!get.bdcm.w.fixed()) 
    { add.general.scaleOperator(bdcCollapseWtScalerID(), bdcCollapseWtID(), 0.75, 1) }
  
  add.comment("Scalers for smcTree model (hyper)parameters") 
  if (!get.smct.popsf.fixed()) 
   { add.general.scaleOperator(popSFScalerID(), popSFID(), 0.75, 1) }
  add.general.scaleOperator(bdcOriginHtScalerID(), bdcOriginHtID(), 0.75, 1)
  
  pergtreewt <- 60 / nof.alignments()
  add.bigcomment("Operators for the trees for each locus")
  gtrees <- get.gtrees()
  for (g in 1:length(gtrees)) {
    a <- get.gtrees()[[g]]$first.partition
    if (length(TheAlignmentData$all.sequences[[a]]) > 2) { 
      add.comment(paste0("topology-changing ops for tree for locus ", g)) 
      add.subtreeSlide(g, gtreewweight(15)) 
      add.narrowExchange(g, gtreewweight(15)) 
      add.wideExchange(g, gtreewweight(4)) 
      add.wilsonBalding(g, gtreewweight(4)) 
    }
    add.comment(paste0("height-only-changing ops for tree for locus ", g))
    add.uniformOperator.internalNodeHeights(g, gtreewweight(30)) 
    add.scaleOperator.allHeights(g, gtreewweight(4)) 
    add.scaleOperator.rootHeight(g, gtreewweight(4)) 
    add.comment(paste0("op for heights and relative clock rate for for tree for locus ", g))
    add.upRate.downHeights(g, gtreewweight(30))
  }
  
  add.bigcomment("Operators for clock rates and subst models")
  clocks <- get.clocks()
  if (length(clocks) > 1) {
    for (c in 2:length(clocks)) {
      add.comment(paste0("op for relative clock rate for partition ", c)) 
      add.general.scaleOperator(clockRateScalerID(c), clockRateID.c(c), 0.75, gtreewweight(2))
    }
  }
  
  siteMs <- get.siteMs() 
  for (u in 1:length(siteMs)) {
    add.comment(paste0("kappas and frequencies ops for subst model ", u))
    #TODO other subst models
    add.general.scaleOperator(kappaScalerID.u(u), kappaID.u(u), 0.5, gtreewweight(2))
    add.deltaExchange.frequencies(u, gtreewweight(2))
  }
}



#################################################################################
# low level, more general functions

add.nodeReheight <- function(wt) {
  attrs <- c(id=nodeReheightID(), spec="NodeReheight", 
             taxonset=IDtoREF(taxonSetOfSetsID()), tree=IDtoREF(smcTreeID()), weight=wt)
  add.opennode("operator", attrs=attrs)
  gtrees <- get.gtrees()
  for (g in 1:length(gtrees)) {
    add.node("tree", attrs=c(idref=geneTreeID.g(g), name="genetree"))
  }
  add.closetag()
} 


add.nodesNudge <- function(wt) {
  attrs <- c(id=nodesNudgeID(), spec="stacey.NodesNudge", weight=wt)
  add.opennode("operator", attrs=attrs)
  add.node("smcTree", attrs=c(idref=smcTreeID()))
  gtrees <- get.gtrees()
  for (g in 1:length(gtrees)) {
    add.node("tree", attrs=c(idref=geneTreeID.g(g), name="geneTree"))
  }
  add.closetag()
}


add.coordinatedPruneRegraft <- function(wt) {
  attrs <- c(id=coordinatedPruneRegraftID(), spec="stacey.CoordinatedPruneRegraft", weight=wt)
  add.opennode("operator", attrs=attrs)
  add.node("smcTree", attrs=c(idref=smcTreeID()))
  gtrees <- get.gtrees()
  for (g in 1:length(gtrees)) {
    add.node("tree", attrs=c(idref=geneTreeID.g(g), name="geneTree"))
  }
  add.closetag()
  
}


add.focusedScaler <- function(wt) {
  attrs <- c(id=focusedScalerID(), spec="stacey.FocusedNodeHeightScaler", weight=wt)
  add.opennode("operator", attrs=attrs)
  add.node("smcTree", attrs=c(idref=smcTreeID()))
  gtrees <- get.gtrees()
  for (g in 1:length(gtrees)) {
    add.node("tree", attrs=c(idref=geneTreeID.g(g), name="geneTree"))
  }
  add.closetag()
}


add.general.scaleOperator <- function(ID, paramID, sfactor, weight) {
  attrs <- c(id=ID, spec="ScaleOperator", parameter=IDtoREF(paramID), scaleFactor=sfactor, weight=weight)
  add.node("operator", attrs=attrs)
}  


gtreewweight <- function(perc) {
  pergtreewt <- 60 / nof.alignments()
  ceiling(100 * perc * pergtreewt) / 100
}


#################################################################################


add.upGrowthClocks.downPopsHeights <- function(wt) {
  gtrees <- get.gtrees()
  clocks <- get.clocks()
  bdcm <- get.bdc.model()
  attrs <- c(id=upGrowthClocks.downPopsHeightsID(), spec="UpDownOperator", weight=wt, scaleFactor="0.75")
  add.opennode("operator", attrs=attrs)
  
  if (!get.bdcm.growthrate.fixed()) {
    add.node("parameter", attrs=c(idref=bdcGrowthID(), name="up"))
  }
  if (length(clocks) > 1) {
    for (c in 2:length(clocks)) {
      add.node("parameter", attrs=c(idref=clockRateID.c(c), name="up"))
    }
  }
  if (!get.smct.popsf.fixed()) {
    add.node("parameter", attrs=c(idref=popSFID(), name="down"))
  }
  add.node("tree", attrs=c(idref=smcTreeID(), name="down"))
  
  for (g in 1:length(gtrees)) {
    add.node("tree", attrs=c(idref=geneTreeID.g(g), name="down"))
  }  
  add.closetag() 
} 

###################### gtree ops ###########################################


add.subtreeSlide <- function(g, weight) {
  add.node("operator", attrs=c(id=subtreeSlideID(g), spec="SubtreeSlide", 
                               tree=IDtoREF(geneTreeID.g(g)), weight=weight))
}


add.narrowExchange <- function(g, weight) {
  add.node("operator", attrs=c(id=narrowID(g), spec="Exchange", 
                               tree=IDtoREF(geneTreeID.g(g)), weight=weight))
}


add.wideExchange <- function(g, weight) {
  add.node("operator", attrs=c(id=wideID(g), spec="Exchange", 
                               tree=IDtoREF(geneTreeID.g(g)), isNarrow="false", weight=weight))
}


add.wilsonBalding <- function(g, weight) {
  add.node("operator", attrs=c(id=WilsonBaldingID(g), spec="WilsonBalding", 
                               tree=IDtoREF(geneTreeID.g(g)), weight=weight))
}


add.uniformOperator.internalNodeHeights <- function(g, weight) {
  add.node("operator", attrs=c(id=uniformID(g), spec="Uniform",
                               tree=IDtoREF(geneTreeID.g(g)), weight=weight))
}


add.scaleOperator.rootHeight <- function(g, weight) {
  add.node("operator", attrs=c(id=treeRootScalerID(g), spec="ScaleOperator", 
                               tree=IDtoREF(geneTreeID.g(g)), scaleFactor="0.5", weight=weight))
}


add.scaleOperator.allHeights <- function(g, weight) {
  add.node("operator", attrs=c(id=treeScalerID(g), spec="ScaleOperator",
                               tree=IDtoREF(geneTreeID.g(g)), scaleFactor="0.5", weight=weight))
}


#TODO I'm confused here. which clocks, which gtrees?
add.upRate.downHeights <- function(a, weight) {
  attrs <- c(id=upClock.downHeightsID(a), spec="UpDownOperator", scaleFactor="0.75", weight=weight)
  add.opennode("operator", attrs=attrs)
  if (a > 1) {
    add.node("parameter", attrs=c(idref=clockRateID.a(a), name="up"))
  } 
  add.node("tree", attrs=c(idref=geneTreeID.a(a), name="down"))
  add.closetag()
}




add.deltaExchange.frequencies <- function(u, weight) {
  attrs <- c(id=frequenciesExchangerID(u), spec="DeltaExchangeOperator", delta="0.01", weight=weight)
  add.opennode("operator", attrs=attrs)
  add.node("parameter", attrs=c(idref=frequenciesParamID.u(u)))
  add.closetag() 
}





