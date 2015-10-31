
add.operators <- function(A) {
  num.mincls <- length(get.minclusters())
  num.gtrees <- length(get.GTrees())
  
  add.comment("Operators for both smcTree and locus trees.")
  stree.op.wt <- signif(400 / sqrt(num.gtrees), digits=3)
  delay.multiplier <- 0
  add.staceyNodeReheight(3*stree.op.wt, delay.multiplier)
  if (num.mincls > 2) {
    add.ThreeBranchAdjuster(stree.op.wt, delay.multiplier)
    add.nodesNudge(stree.op.wt, delay.multiplier)
    add.coordinatedPruneRegraft(stree.op.wt, delay.multiplier)
  }
  if (num.mincls > 4) {
    add.focusedScaler(stree.op.wt, delay.multiplier)
  }
  
  add.comment("stretch and squeeze everything") 
  add.upGrowthClocks.downPopsHeights(2, "0p100", 0.1) 
  add.upGrowthClocks.downPopsHeights(2, "0p500", 0.5) 
  add.upGrowthClocks.downPopsHeights(2, "0p850", 0.85) 
  add.upGrowthClocks.downPopsHeights(2, "0p900", 0.9) 
  add.upGrowthClocks.downPopsHeights(2, "0p970", 0.97) 
  add.upGrowthClocks.downPopsHeights(2, "0p990", 0.99) 
  add.upGrowthClocks.downPopsHeights(2, "0p997", 0.997) 
  add.upGrowthClocks.downPopsHeights(2, "0p999", 0.999) 
  
  add.comment("Scalers for the birth-death-collapse model (hyper)parameters") 
  
  bdc <- get.bdc.model()
  bdc.g <- get.child(bdc, "growthrate")
  bdc.rd <- get.child(bdc, "reldeath")
  bdc.w <- get.child(bdc, "w")
  if (is.estimated(bdc.g)) {
    add.general.no.opt.scaleOperator(bdcGrowthRateScalerID("0p100"), bdcGrowthID(), 0.1, .4)
    add.general.no.opt.scaleOperator(bdcGrowthRateScalerID("0p500"), bdcGrowthID(), 0.5, .4)
    add.general.no.opt.scaleOperator(bdcGrowthRateScalerID("0p900"), bdcGrowthID(), 0.9, .4)
    add.general.no.opt.scaleOperator(bdcGrowthRateScalerID("0p990"), bdcGrowthID(), 0.99, .4)
    add.general.no.opt.scaleOperator(bdcGrowthRateScalerID("0p999"), bdcGrowthID(), 0.999, .4)
  }
  if (is.estimated(bdc.rd)) {
    add.general.no.opt.scaleOperator(bdcRelDeathScalerID("0p010"), bdcRelDeathID(), 0.01, .2)
    add.general.no.opt.scaleOperator(bdcRelDeathScalerID("0p100"), bdcRelDeathID(), 0.1, .2)
    add.general.no.opt.scaleOperator(bdcRelDeathScalerID("0p500"), bdcRelDeathID(), 0.5, .2)
    add.general.no.opt.scaleOperator(bdcRelDeathScalerID("0p900"), bdcRelDeathID(), 0.9, .2)
    add.general.no.opt.scaleOperator(bdcRelDeathScalerID("0p990"), bdcRelDeathID(), 0.99, .2)
    
    add.general.runifrandwalkOperator(bdcRelDeathRandWalkerID("0p0001"), bdcRelDeathID(), 0.0001,   .2)
    add.general.runifrandwalkOperator(bdcRelDeathRandWalkerID("0p0010"), bdcRelDeathID(), 0.001,    .2)
    add.general.runifrandwalkOperator(bdcRelDeathRandWalkerID("0p0100"), bdcRelDeathID(), 0.01,     .2)
    add.general.runifrandwalkOperator(bdcRelDeathRandWalkerID("0p1000"), bdcRelDeathID(), 0.1,      .2)
    add.general.runifrandwalkOperator(bdcRelDeathRandWalkerID("0p3000"), bdcRelDeathID(), 0.3,      .2)
  }
  if (is.estimated(bdc.w)) {
    add.general.no.opt.scaleOperator(bdcCollapseWtScalerID("0p0001"), bdcCollapseWtID(), 0.0001, .2)
    add.general.no.opt.scaleOperator(bdcCollapseWtScalerID("0p0010"), bdcCollapseWtID(), 0.001, .2)
    add.general.no.opt.scaleOperator(bdcCollapseWtScalerID("0p0100"), bdcCollapseWtID(), 0.01, .2)
    add.general.no.opt.scaleOperator(bdcCollapseWtScalerID("0p1000"), bdcCollapseWtID(), 0.1, .2)
    add.general.no.opt.scaleOperator(bdcCollapseWtScalerID("0p5000"), bdcCollapseWtID(), 0.5, .2)
    add.general.no.opt.scaleOperator(bdcCollapseWtScalerID("0p9000"), bdcCollapseWtID(), 0.9, .2)
    
    add.general.runifrandwalkOperator(bdcCollapseWtRandWalkerID("0p000001"), bdcCollapseWtID(), 0.000001, .2)
    add.general.runifrandwalkOperator(bdcCollapseWtRandWalkerID("0p00001"), bdcCollapseWtID(), 0.00001,  .2)
    add.general.runifrandwalkOperator(bdcCollapseWtRandWalkerID("0p0001"), bdcCollapseWtID(), 0.0001,   .2)
    add.general.runifrandwalkOperator(bdcCollapseWtRandWalkerID("0p001"), bdcCollapseWtID(), 0.001,    .2)
    add.general.runifrandwalkOperator(bdcCollapseWtRandWalkerID("0p010"), bdcCollapseWtID(), 0.01,     .2)
    add.general.runifrandwalkOperator(bdcCollapseWtRandWalkerID("0p100"), bdcCollapseWtID(), 0.1,      .2)
  }
  
  add.comment("Scalers for smcTree model (hyper)parameters") 
  
  popSF <- get.child(get.smc.coalescent(), "popSF") 
  if (is.estimated(popSF)) {
    add.general.no.opt.scaleOperator(popSFScalerID("0p1"), popSFID(), 0.1, 1)
    add.general.no.opt.scaleOperator(popSFScalerID("0p5"), popSFID(), 0.5, 1)
    add.general.no.opt.scaleOperator(popSFScalerID("0p9"), popSFID(), 0.9, 1)
    add.general.no.opt.scaleOperator(popSFScalerID("0p99"), popSFID(), 0.99, 1)
    add.general.no.opt.scaleOperator(popSFScalerID("0p999"), popSFID(), 0.999, 1)
  }
  
  add.general.no.opt.scaleOperator(bdcOriginHeightScalerID("0p1"), bdcOriginHtID(), 0.1, .4)
  add.general.no.opt.scaleOperator(bdcOriginHeightScalerID("0p5"), bdcOriginHtID(), 0.5, .4)
  add.general.no.opt.scaleOperator(bdcOriginHeightScalerID("0p9"), bdcOriginHtID(), 0.9, .4)
  add.general.no.opt.scaleOperator(bdcOriginHeightScalerID("0p99"), bdcOriginHtID(), 0.99, .4)
  add.general.no.opt.scaleOperator(bdcOriginHeightScalerID("0p999"), bdcOriginHtID(), 0.999, .4)
  
  pergtreewt <- 60 / get.nof.alignments()
  add.bigcomment("Operators for the trees for each locus")
  for (g in 1:num.gtrees) {
    a <- get.GTrees.with.first.partition()[[g]]$first.partition
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
  clocks <- get.Clocks()
  if (length(clocks) > 1) {
    for (u in 2:length(clocks)) {
    	prm <- get.child(clocks[[u]], "PartitionRateM")
    	clkrate <- get.child(prm, "rate")
      if (is.estimated(clkrate)) {
        add.comment(paste0("op for relative clock rate for partition ", u)) 
        add.general.scaleOperator(clockRateScalerID(u), clockRateID.u(u), 0.75, gtreewweight(2))        
      }
     }
  }
  
  siteMs <- get.SiteMs() 
  for (u in 1:length(siteMs)) {
    add.comment(paste0("kappas and frequencies ops for subst model ", u))
    #TODO other subst models
    substM <- get.child(siteMs[[u]], "SubstM")
    model <- get.child(substM, "model")
    if (model$kind == "HKY") {
    	kappa <- get.childvalue(model, "kappa")
    	freqs <- get.childvalue(model, "freqs")
    	if (is.numeric(kappa)) {
        # no operator
      } else {
        add.general.scaleOperator(kappaScalerID.u(u), kappaID.u(u), 0.5, gtreewweight(2))
      }      
      if (is.character(freqs)) {
        # no operator
      } else {
        add.deltaExchange.frequencies(u, gtreewweight(2))
      }
    }
  }
}



#################################################################################
# low level, more general functions

add.staceyNodeReheight <- function(wt, delay.multiplier) {
  dly <- as.character(as.integer(wt*delay.multiplier))
  attrs <- c(id=nodeReheightID(), spec="stacey.StaceyNodeReheight", 
             proportionUniform="0.05", delay=dly, weight=wt)
  open.xmlnode("operator", attrs=attrs)
  add.xmlnode("popSF", attrs=c(idref=popSFID()))
  add.xmlnode("smcTree", attrs=c(idref=smcTreeID()))
  gtrees <- get.GTrees()
  for (u in 1:length(gtrees)) {
    add.xmlnode("tree", attrs=c(idref=geneTreeID.u(u), name="geneTree"))
  }
  close.xmlnode()
} 



add.ThreeBranchAdjuster <- function(wt, delay.multiplier) {
  dly <- as.character(as.integer(wt*delay.multiplier))
  attrs <- c(id=threeBranchAdjusterID(), spec="stacey.ThreeBranchAdjuster", delay=dly, weight=wt)
  open.xmlnode("operator", attrs=attrs)
  add.xmlnode("popSF", attrs=c(idref=popSFID()))
  add.xmlnode("smcTree", attrs=c(idref=smcTreeID()))
  gtrees <- get.GTrees()
  for (u in 1:length(gtrees)) {
    add.xmlnode("tree", attrs=c(idref=geneTreeID.u(u), name="geneTree"))
  }
  close.xmlnode()
}



add.nodesNudge <- function(wt, delay.multiplier) {
  dly <- as.character(as.integer(wt*delay.multiplier))
  attrs <- c(id=nodesNudgeID(), spec="stacey.NodesNudge", delay=dly, weight=wt)
  open.xmlnode("operator", attrs=attrs)
  add.xmlnode("smcTree", attrs=c(idref=smcTreeID()))
  gtrees <- get.GTrees()
  for (u in 1:length(gtrees)) {
    add.xmlnode("tree", attrs=c(idref=geneTreeID.u(u), name="geneTree"))
  }
  close.xmlnode()
}


add.coordinatedPruneRegraft <- function(wt, delay.multiplier) {
  dly <- as.character(as.integer(wt*delay.multiplier))
  attrs <- c(id=coordinatedPruneRegraftID(), spec="stacey.CoordinatedPruneRegraft", delay=dly, weight=wt)
  open.xmlnode("operator", attrs=attrs)
  add.xmlnode("smcTree", attrs=c(idref=smcTreeID()))
  gtrees <- get.GTrees()
  for (u in 1:length(gtrees)) {
    add.xmlnode("tree", attrs=c(idref=geneTreeID.u(u), name="geneTree"))
  }
  close.xmlnode()
  
}


add.focusedScaler <- function(wt, delay.multiplier) {
  dly <- as.character(as.integer(wt*delay.multiplier))
  attrs <- c(id=focusedScalerID(), spec="stacey.FocusedNodeHeightScaler", delay=dly, weight=wt)
  open.xmlnode("operator", attrs=attrs)
  add.xmlnode("smcTree", attrs=c(idref=smcTreeID()))
  gtrees <- get.GTrees()
  for (u in 1:length(gtrees)) {
    add.xmlnode("tree", attrs=c(idref=geneTreeID.u(u), name="geneTree"))
  }
  close.xmlnode()
}


# 
# add.heightsWarper <- function(wt) {
#   attrs <- c(id=heightsWarperID(), spec="stacey.HeightsWarper", weight=wt)
#   open.xmlnode("operator", attrs=attrs)
#   add.xmlnode("smcTree", attrs=c(idref=smcTreeID()))
#   gtrees <- get.gtrees()
#   for (u in 1:length(gtrees)) {
#     add.xmlnode("geneTree", attrs=c(idref=geneTreeCoalFactorID.u(u)))
#   }
#   close.xmlnode()
#   
# }




add.general.scaleOperator <- function(ID, paramID, sfactor, weight) {
  attrs <- c(id=ID, spec="ScaleOperator", parameter=IDtoREF(paramID), scaleFactor=sfactor, weight=weight)
  add.xmlnode("operator", attrs=attrs)
}  

add.general.no.opt.scaleOperator <- function(ID, paramID, sfactor, weight) {
  attrs <- c(id=ID, spec="ScaleOperator", parameter=IDtoREF(paramID), scaleFactor=sfactor, weight=weight, optimise="false")
  add.xmlnode("operator", attrs=attrs)
}  


add.general.runifrandwalkOperator <- function(ID, paramID, halfwidth, weight) {
  attrs <- c(id=ID, spec="RealRandomWalkOperator", parameter=IDtoREF(paramID), windowSize=halfwidth, weight=weight)
  add.xmlnode("operator", attrs=attrs)
}  


gtreewweight <- function(perc) {
  pergtreewt <- 60 / get.nof.alignments()
  ceiling(100 * perc * pergtreewt) / 100
}


#################################################################################


add.upGrowthClocks.downPopsHeights <- function(wt, sf.text, sf) {
  gtrees <- get.GTrees()
  clocks <- get.Clocks()
  bdc <- get.bdc.model()
  bdc.g <- get.child(bdc, "growthrate")
  attrs <- c(id=upGrowthClocks.downPopsHeightsID(sf.text), spec="UpDownOperator", weight=wt, scaleFactor=sf, optimise="false")
  open.xmlnode("operator", attrs=attrs)
  
  if (is.estimated(bdc.g)) {
    add.xmlnode("parameter", attrs=c(idref=bdcGrowthID(), name="up"))
  }

  for (u in 1:length(clocks)) {
  	prm <- get.child(clocks[[u]], "PartitionRateM")
  	clkrate <- get.child(prm, "rate")
    if (is.estimated(clkrate)) {
      add.xmlnode("parameter", attrs=c(idref=clockRateID.u(u), name="up"))
    }
  }

  popSF <- get.child(get.smc.coalescent(), "popSF")
  if (is.estimated(popSF)) {
    add.xmlnode("parameter", attrs=c(idref=popSFID(), name="down"))
  }
  add.xmlnode("tree", attrs=c(idref=smcTreeID(), name="down"))
  
  for (u in 1:length(gtrees)) {
    add.xmlnode("tree", attrs=c(idref=geneTreeID.u(u), name="down"))
  }  
  close.xmlnode() 
} 

###################### gtree ops ###########################################


add.subtreeSlide <- function(u, weight) {
  add.xmlnode("operator", attrs=c(id=subtreeSlideID(u), spec="SubtreeSlide", 
                               tree=IDtoREF(geneTreeID.u(u)), weight=weight))
}


add.narrowExchange <- function(u, weight) {
  add.xmlnode("operator", attrs=c(id=narrowID(u), spec="Exchange", 
                               tree=IDtoREF(geneTreeID.u(u)), weight=weight))
}


add.wideExchange <- function(u, weight) {
  add.xmlnode("operator", attrs=c(id=wideID(u), spec="Exchange", 
                               tree=IDtoREF(geneTreeID.u(u)), isNarrow="false", weight=weight))
}


add.wilsonBalding <- function(u, weight) {
  add.xmlnode("operator", attrs=c(id=WilsonBaldingID(u), spec="WilsonBalding", 
                               tree=IDtoREF(geneTreeID.u(u)), weight=weight))
}


add.uniformOperator.internalNodeHeights <- function(u, weight) {
  add.xmlnode("operator", attrs=c(id=uniformID(u), spec="Uniform",
                               tree=IDtoREF(geneTreeID.u(u)), weight=weight))
}


add.scaleOperator.rootHeight <- function(u, weight) {
  add.xmlnode("operator", attrs=c(id=treeRootScalerID(u), spec="ScaleOperator", 
                               tree=IDtoREF(geneTreeID.u(u)), scaleFactor="0.5", weight=weight))
}


add.scaleOperator.allHeights <- function(u, weight) {
  add.xmlnode("operator", attrs=c(id=treeScalerID(u), spec="ScaleOperator",
                               tree=IDtoREF(geneTreeID.u(u)), scaleFactor="0.5", weight=weight))
}


#TODO I'm confused here. which clocks, which gtrees?
add.upRate.downHeights <- function(a, weight) {
  attrs <- c(id=upClock.downHeightsID(a), spec="UpDownOperator", scaleFactor="0.75", weight=weight)
  open.xmlnode("operator", attrs=attrs)
  if (a > 1) {
    add.xmlnode("parameter", attrs=c(idref=clockRateID.a(a), name="up"))
  } 
  add.xmlnode("tree", attrs=c(idref=geneTreeID.a(a), name="down"))
  close.xmlnode()
}




add.deltaExchange.frequencies <- function(u, weight) {
  attrs <- c(id=frequenciesExchangerID(u), spec="DeltaExchangeOperator", delta="0.01", weight=weight)
  open.xmlnode("operator", attrs=attrs)
  add.xmlnode("parameter", attrs=c(idref=frequenciesParamID.u(u)))
  close.xmlnode() 
}





