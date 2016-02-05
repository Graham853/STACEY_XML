
add.operators <- function() {
  num.mincls <- length(get.minclusters())
  num.gtrees <- length(get.GTrees())
  
  wt.Reheight <- as.numeric(get.runoption("op.wt.Reheight"))
  wt.Nudge <-    as.numeric(get.runoption("op.wt.Nudge"))
  wt.Focused <-  as.numeric(get.runoption("op.wt.Focused"))
  wt.Regraft <-  as.numeric(get.runoption("op.wt.Regraft"))
  wt.BAdjust <-  as.numeric(get.runoption("op.wt.BAdjust"))
  
  
  add.comment("Operators for both smcTree and locus trees.")
  stree.op.wt <- signif(400 / sqrt(num.gtrees), digits=3)
  delay.multiplier <- 0
  add.staceyNodeReheight(3*stree.op.wt*wt.Reheight, delay.multiplier)
  if (num.mincls > 2) {
    add.ThreeBranchAdjuster(stree.op.wt*wt.BAdjust, delay.multiplier)
    add.nodesNudge(stree.op.wt*wt.Nudge, delay.multiplier)
    add.coordinatedPruneRegraft(stree.op.wt*wt.Regraft, delay.multiplier)
  }
  if (num.mincls > 4) {
    add.focusedScaler(stree.op.wt*wt.Focused, delay.multiplier)
  }
  
  add.comment("stretch and squeeze everything") 
  add.upGrowthClocks.downPopsHeights(2, "0p100", "0.100") 
  add.upGrowthClocks.downPopsHeights(2, "0p350", "0.350") 
  add.upGrowthClocks.downPopsHeights(2, "0p700", "0.700") 
  add.upGrowthClocks.downPopsHeights(2, "0p900", "0.900") 
  add.upGrowthClocks.downPopsHeights(2, "0p970", "0.970") 
  add.upGrowthClocks.downPopsHeights(2, "0p990", "0.990") 
  add.upGrowthClocks.downPopsHeights(2, "0p997", "0.997") 
  add.upGrowthClocks.downPopsHeights(2, "0p999", "0.999") 
  
  
  add.comment("Scalers for the birth-death-collapse model (hyper)parameters") 
  
  bdc <- get.bdc.model()
  bdc.g <- get.child(bdc, "growthrate")
  bdc.rd <- get.child(bdc, "reldeath")
  bdc.w <- get.child(bdc, "w")
  scales <- c("100", "350", "700", "900", "970", "990", "997", "999")
  halfwidths <- c("001", "010", "030", "100", "300")
  if (is.estimated(bdc.g)) {
    add.set.of.scaleOperators(scales, bdcGrowthRateScalerIDroot(), bdcGrowthID(), 0.3)
    }
  if (is.estimated(bdc.rd)) {
    add.set.of.scaleOperators(scales, bdcRelDeathScalerIDroot(), bdcRelDeathID(), 0.1)
    add.set.of.runifrandwalkOperators(halfwidths, bdcRelDeathRandWalkerIDroot(), bdcRelDeathID(), 0.2)
  }
  if (is.estimated(bdc.w)) {
    add.set.of.scaleOperators(scales, bdcCollapseWtScalerIDroot(), bdcCollapseWtID(), 0.1)
    add.set.of.runifrandwalkOperators(halfwidths, bdcCollapseWtRandWalkerIDroot(), bdcCollapseWtID(), 0.2)
  }
  
  add.comment("Scalers for smcTree model (hyper)parameters") 
  
  popSF <- get.child(get.smc.coalescent(), "popSF") 
  if (is.estimated(popSF)) {
    add.set.of.scaleOperators(scales, popSFScalerIDroot(), popSFID(), 0.7)
  }
  add.set.of.scaleOperators(scales, bdcOriginHeightScalerIDroot(), bdcOriginHtID(), 0.3)
  
  add.bigcomment("Operators for the trees for each locus")
  for (g in 1:num.gtrees) {
    a <- get.GTrees.with.first.partition()[[g]]$first.partition
    if (length(TheAlignmentData$all.sequences[[a]]) > 2) { 
      add.comment(paste0("topology-changing ops for tree for locus ", g)) 
      add.setof.subtreeSlide(g, gtreewweight(30)) 
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
  if (wt > 0) {
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
} 



add.ThreeBranchAdjuster <- function(wt, delay.multiplier) {
  if (wt > 0) {
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
}



add.nodesNudge <- function(wt, delay.multiplier) {
  if (wt > 0) {
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
}


add.coordinatedPruneRegraft <- function(wt, delay.multiplier) {
  if (wt > 0) {
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
}



add.focusedScaler <- function(wt, delay.multiplier) {
  if (wt > 0) {
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
}



add.general.scaleOperator <- function(ID, paramID, sfactor, weight) {
  attrs <- c(id=ID, spec="ScaleOperator", parameter=IDtoREF(paramID), scaleFactor=sfactor, weight=weight)
  add.xmlnode("operator", attrs=attrs)
}  

add.general.no.opt.scaleOperator <- function(ID, paramID, sfactor, weight) {
  attrs <- c(id=ID, spec="ScaleOperator", parameter=IDtoREF(paramID), scaleFactor=sfactor, weight=weight, optimise="false")
  add.xmlnode("operator", attrs=attrs)
}  

add.set.of.scaleOperators <- function(scales, OpIDroot, paramID, weight) {
  for (i in 1:length(scales)) {
    opID <- paste0(OpIDroot, scales[i])
    scale <- sprintf("%1.3f", as.double(scales[i])/1000)
    add.general.no.opt.scaleOperator(opID, paramID, scale, weight)
  }
}


add.general.runifrandwalkOperator <- function(ID, paramID, halfwidth, weight) {
  attrs <- c(id=ID, spec="RealRandomWalkOperator", parameter=IDtoREF(paramID), windowSize=halfwidth, weight=weight)
  add.xmlnode("operator", attrs=attrs)
}  


add.set.of.runifrandwalkOperators <- function(halfwidths, OpIDroot, paramID, weight) {
  for (i in 1:length(halfwidths)) {
    opID <- paste0(OpIDroot, halfwidths[i])
    halfwidth <- sprintf("%1.3f", as.double(halfwidths[i])/1000)
    add.general.runifrandwalkOperator(opID, paramID, halfwidth, weight)
  }
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


add.setof.subtreeSlide <- function(u, weight) {
  sizes <- c("0001", "0003", "0010", "0030", "0100", "0300", "1000", "3000")
  for (i in 1:length(sizes)) {
    opID <- paste0(subtreeSlideIDroot(u), sizes[i])
    size <- sprintf("%1.4f", as.double(sizes[i])/10000)
    attrs <- c(id=opID, spec="SubtreeSlide", tree=IDtoREF(geneTreeID.u(u)),
               size=size, optimise="false", weight=weight/2/length(sizes))
    add.xmlnode("operator", attrs=attrs)
  }
  opID <- paste0(subtreeSlideIDroot(u), "default")
  attrs <- c(id=opID, spec="SubtreeSlide", tree=IDtoREF(geneTreeID.u(u)),
             weight=weight)
    add.xmlnode("operator", attrs=attrs)
}
#public Input<Double> sizeInput = new Input<Double>("size", "size of the slide, default 1.0", 1.0);
#public Input<Boolean> gaussianInput = new Input<Boolean>("gaussian", "Gaussian (=true=default) or uniform delta", true);
#public Input<Boolean> optimiseInput = new Input<Boolean>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);




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





