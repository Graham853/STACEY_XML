
add.operators <- function(A) {
  nofmcs <- length(get.minclusters())
  add.comment("Operators for both smcTree and gene trees.")
  add.smctree.gtree.op(18, nodeReheightID(),            "NodeReheight")
  if (nofmcs > 2)
  { add.smctree.gtree.op(6,  nodesNudgeID(),              "stacey.NodesNudge") }
  if (nofmcs > 2)
  { add.smctree.gtree.op(6,  coordinatedPruneRegraftID(), "stacey.CoordinatedPruneRegraft") }
  if (nofmcs > 3)
  { add.smctree.gtree.op(6,  focusedScalerID(),           "stacey.FocusedNodeHeightScaler") }
  
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
  add.bigcomment("Operators for the gene trees")
  for (g in 1:nof.alignments()) {
    if (nofmcs > 2) { 
      #TODO should be nof sequences for gene
      add.comment(paste0("topology-changing ops for gene tree ", g)) 
      add.subtreeSlide(g, gtreewweight(15)) 
      add.narrowExchange(g, gtreewweight(15)) 
      add.wideExchange(g, gtreewweight(4)) 
      add.wilsonBalding(g, gtreewweight(4)) 
    }
    add.comment(paste0("height-only-changing ops for gene tree ", g))
    add.uniformOperator.internalNodeHeights(g, gtreewweight(30)) 
    add.scaleOperator.allHeights(g, gtreewweight(4)) 
    add.scaleOperator.rootHeight(g, gtreewweight(4)) 
    add.comment(paste0("op for heights and relative clock rate for gene tree ", g))
    add.upRate.downHeights(g, gtreewweight(30))
  }
  
  add.bigcomment("Operators for clock rates and subst models")
  clocks <- get.clocks()
  if (length(clocks) > 1) {
    for (c in 2:length(clocks)) {
      add.comment(paste0("op for relative clock rate for gene tree ", c)) 
      add.general.scaleOperator(clockRateScalerID(c), clockRateID(c), 0.75, gtreewweight(2))
    }
  }
  
  substs <- get.substs() 
  for (s in 1:length(substs)) {
    add.comment(paste0("kappas and frequencies ops for subst model ", s))
    #TODO other subst models
    add.general.scaleOperator(kappaScalerID(s), kappaID(s), 0.5, gtreewweight(2))
    add.deltaExchange.frequencies(s, gtreewweight(2))
  }
  
  



}



#################################################################################
# low level, more general functions

add.smctree.gtree.op <- function(wt, id, spec) {
  attrs <- c(id=id, spec=spec, taxonset=IDtoREF(taxonSetOfSetsID()), tree=IDtoREF(smcTreeID()), weight=wt)
  add.opennode("operator", attrs=attrs)
  for (g in 1:nof.alignments()) {
    add.node("tree", c(id=geneTreeID(g), name="genetree"))
  }
  add.closetag()
} 



add.general.scaleOperator <- function(ID, param, sfactor, weight) {
  attrs <- c(id=ID, spec="ScaleOperator", scaleFactor=sfactor, weight=weight)
  add.opennode("operator", attrs=attrs)
  add.node("parameter", id=param)
  add.closetag()
}  


gtreewweight <- function(perc) {
  pergtreewt <- 60 / nof.alignments()
  ceiling(100 * perc * pergtreewt) / 100
}


#################################################################################


add.upGrowthClocks.downPopsHeights <- function(wt) {
  bdcm <- get.bdc.model()
  attrs <- c(id=upGrowthClocks.downPopsHeightsID(), spec="UpDownOperator", weight=wt, scaleFactor="0.75")
  add.opennode("operator", attrs=attrs)
  
  if (!get.bdcm.growthrate.fixed()) {
    add.node("parameter", attrs=c(id=bdcGrowthID(), name="up"))
  }
  if (nof.alignments() > 1) {
    for (g in 2:nof.alignments()) {
      add.node("parameter", attrs=c(id=clockRateID(g), name="up"))
    }
  }
  if (!get.smct.popsf.fixed()) {
    add.node("parameter", attrs=c(id=popSFID(), name="down"))
  }
  add.node("tree", attrs=c(id=smcTreeID(), name="down"))
  
  for (g in 1:nof.alignments()) {
    add.node("tree", attrs=c(id=geneTreeID(g), name="down"))
  }  
  add.closetag() 
} 

###################### gtree ops ###########################################


add.subtreeSlide <- function(g, weight) {
  add.node("operator", attrs=c(id=subtreeSlideID(g), spec="SubtreeSlide", 
                               tree=IDtoREF(geneTreeID(g)), weight=weight))
}


add.narrowExchange <- function(g, weight) {
  add.node("operator", attrs=c(id=narrowID(g), spec="Exchange", 
                               tree=IDtoREF(geneTreeID(g)), weight=weight))
}


add.wideExchange <- function(g, weight) {
  add.node("operator", attrs=c(id=wideID(g), spec="Exchange", 
                               tree=IDtoREF(geneTreeID(g)), isNarrow="false", weight=weight))
}


add.wilsonBalding <- function(g, weight) {
  add.node("operator", attrs=c(id=WilsonBaldingID(g), spec="WilsonBalding", 
                               tree=IDtoREF(geneTreeID(g)), weight=weight))
}


add.uniformOperator.internalNodeHeights <- function(g, weight) {
  add.node("operator", attrs=c(id=uniformID(g), spec="Uniform",
                               tree=IDtoREF(geneTreeID(g)), weight=weight))
}


add.scaleOperator.rootHeight <- function(g, weight) {
  add.node("operator", attrs=c(id=treeRootScalerID(g), spec="ScaleOperator", 
                               tree=IDtoREF(geneTreeID(g)), scaleFactor="0.5", weight=weight))
}


add.scaleOperator.allHeights <- function(g, weight) {
  add.node("operator", attrs=c(id=treeScalerID(g), spec="ScaleOperator",
                               tree=IDtoREF(geneTreeID(g)), scaleFactor="0.5", weight=weight))
}


add.upRate.downHeights <- function(g, weight) {
  attrs <- c(id=upClock.downHeightsID(g), spec="UpDownOperator", scaleFactor="0.75", weight=weight)
  add.opennode("operator", attrs=attrs)
  if (g > 1) {
    add.node("parameter", attrs=c(id=clockRateID(g), name="up"))
  } 
  add.node("tree", attrs=c(id=geneTreeID(g), name="down"))
  add.closetag()
}




add.deltaExchange.frequencies <- function(g, weight) {
  attrs <- c(id=frequenciesExchangerID(g), spec="DeltaExchangeOperator", delta="0.01", weight=weight)
  add.opennode("operator", attrs=attrs)
  add.node("parameter", attrs=c(id=frequenciesID(g)))
  add.closetag() 
}





