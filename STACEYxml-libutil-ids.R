
######## HELPER FUNCTIONS


IDtoREF <- function(ID) {
  paste0("@", ID)
}



xmlIDFromElement  <- function(x) {
  stopifnot(is.character(x$kind))
  stopifnot(is.character(x$id))
  paste0(x$kind, ".", x$id)
}


xmlIDfromKindID  <- function(kind, id) {
  if (!is.character(kind)) {
    browser()
  }
  
  stopifnot(is.character(kind))
  stopifnot(is.character(id))
  paste0(kind, ".", id)
}


dataID.a <- function(a) {
  fname <- get.gtree.DataFileName(a)$value
  parts <- str_split(fname, fixed("."))[[1]]
  nparts <- length(parts)
  if (nparts == 1) {
    id <- parts[1]
  } else {
    id <- paste(parts[1:(nparts-1)], collapse=".")
  }
  id
}



############ ID FUNCTIONS

partitionID.a <- function(a) {
  paste0("Partition.", id)
}

# data for partition a
partitiondataID.a <- function(a) {
  paste0("Data.", dataID.a(a))
}

# data for first partition belonging gene tree u
partitiondataID.u <- function(u) {
  a <- get.GTrees.with.first.partition()[[u]]$first.partition
  paste0("Data.", dataID.a(a))
}



sequenceID.a <- function(a, taxonname) {
  paste0("Seq.", dataID.a(a), ".", taxonname)
}





initialsmctreeID <- function() {
  "init.smcTree"
}


initialsmctreepopmodelID <- function() {
  "init.smcTree.popmodel"
}


initialsmctreepopsizeID <- function() {
  "init.smcTree.popsize"
}


initialgtreeID <- function(g) {
  paste0("init.gTree", g)
}


initialgtreepopmodelID <- function(g) {
  paste0("init.gTree.popmodel", g)
}


initialgtreepopsizeID <- function(g) {
  paste0("init.gTree.popsize", g)
}


posteriorID <- function() { 
"posterior"
}

likelihoodID <- function() { 
  "likelihood"
}

priorID <- function() { 
  "prior"
}


##########


smcTreeID <- function() {
  smct <- get.child(get.gtree.GTreePrior(1), "prior")
  xmlIDFromElement(smct)
}

taxonSetOfSetsID <- function() { "taxonSetOfSets" }
# spp/min-cluster and sequences IDs not done here

bdcModelID <- function() {
	bdc <- get.bdc.model()
  xmlIDFromElement(bdc)
}


bdcGrowthID <- function() { 
	bdc <- get.bdc.model()
  xmlIDfromKindID(bdc$kind, get.child(bdc, "growthrate")$id)
}


bdcRelDeathID <- function() {
	bdc <- get.bdc.model()
	xmlIDfromKindID(bdc$kind, get.child(bdc, "reldeath")$id)
}


bdcCollapseWtID <- function() {
	bdc <- get.bdc.model()
  xmlIDfromKindID(bdc$kind, get.child(bdc, "w")$id)
}


bdcOriginHtID <- function() {
	bdc <- get.bdc.model()
	xmlIDfromKindID(bdc$kind, get.child(bdc, "oh")$id)
}



smcCoalescentID <- function() {
  xmlIDFromElement(get.smc.coalescent())
}


InvGammaComponentID <- function(c) { 
	smctcoal <- get.smc.coalescent()
  paste0(smctcoal$kind, ".", get.child(smctcoal, "invgammamix")$id, ".", c) 
}


popSFID <- function() {
	smctcoal <- get.smc.coalescent()
  xmlIDfromKindID(smctcoal$kind, get.child(smctcoal, "popSF")$id) 
}



######################################################


########################################################

# GENE TREES. 




geneTreeID.a <- function(a) {
  gtree <- get.gtree.GTree(a)
  xmlIDFromElement(gtree)
}

geneTreeID.u <- function(u) {
  gtree <- get.GTrees()[[u]]
  xmlIDFromElement(gtree)
}



geneTreeBranchRM <- function(u) {
  gtree <- get.GTrees()[[u]]
  xmlIDFromElement(gtree$branchRM$model)
}



geneTreeLhoodID <- function(u) {
  gtree <- get.GTrees()[[u]]
  paste0("GTreeLhood.", gtree$id)
}


geneTreeCoalFactorID.u <- function(u) {
  gtree <- get.GTrees()[[u]]
  paste0("GTreePloidy.", gtree$id)
}


# referred to in init to make random trees
geneTaxonSetID.u <- function(u) {
  gtree <- get.GTrees()[[u]]
  paste0("GTreeTaxonSet.", gtree$id)
}



####################################################
##### SITE MODELS


sitemodelID.a  <- function(a) {
  siteM <- get.gtree.SiteM(a)
  xmlIDFromElement(siteM)
}

sitemodelID.u <- function(u) {
  siteMs <- get.siteMs()
  stopifnot(u <= length(siteMs))
  xmlIDFromElement(siteMs[[u]])
}

# SITE RATE HETEROGENEITY

sitehetMuID.a <- function(a) {
  sitehet <- get.gtree.SiteHet(a)
  stopifnot(is.character(sitehet$value))
  if (sitehet$value == "None") {
    paste0(sitehet$kind, ".mu.", sitehet$id)
  } else {
    #TODO 
  }  
}


sitehetGammaShapeID.a <- function(a) {
  sitehet <- get.gtree.SiteHet(a)
  stopifnot(is.character(sitehet$value))
  if (sitehet$value == "None") {
    paste0(sitehet$kind, ".GammaShape.", sitehet$id)
  } else {
    #TODO 
  }  
}

sitehetPropInvID.a <- function(a) {
  sitehet <- get.gtree.SiteHet(a)
  stopifnot(is.character(sitehet$value))
  if (sitehet$value == "None") {
    paste0(sitehet$kind, ".PropInv.", sitehet$id)
  } else {
    #TODO 
  }  
}




# SUBSTITUTION MODELS 

substsmodelID.a <- function(a) {
  model <- get.child(get.gtree.SubstM(a), "model")
  xmlIDFromElement(model)
}

substsmodelID.u <- function(u) {
	SubstMs <- get.gtree.SubstMs()
	stopifnot(u <= length(SubstMs))
	model <- get.child(SubstMs[[u]], "model")
  xmlIDFromElement(model)
}


kappaID.a <- function(a) {
  model <- get.child(get.gtree.SubstM(a), "model")
  xmlIDfromKindID(model$kind, get.child(model, "kappa")$id)
}

kappaID.u <- function(u) { 
	SubstMs <- get.gtree.SubstMs()
	stopifnot(u <= length(SubstMs))
	model <- get.child(SubstMs[[u]], "model")
	xmlIDfromKindID(model$kind, get.child(model, "kappa")$id)
}


# there are two things that are 'frequencies'. The model and the (vector) parameter.

frequenciesModelID.a <- function(a) {
	model <- get.child(get.gtree.SubstM(a), "model")
	freqs <- get.child(model, "freqs")
  paste0("FrequenciesModel.", freqs$id)
}


frequenciesParamID.a <- function(a) {
  model <- get.child(get.gtree.SubstM(a), "model")
  freqs <- get.child(model, "freqs")
  xmlIDfromKindID(model$kind, freqs$id)
}

frequenciesParamID.u <- function(u) {
	SubstMs <- get.gtree.SubstMs()
	stopifnot(u <= length(SubstMs))
  model <- get.child(SubstMs[[u]], "model")
  xmlIDfromKindID(model$kind, get.child(model, "freqs")$id)
}

empiricalFreqsID.a <- function(a) {
  paste0("empiricalFreqs.", dataID.a(a))
}

# CLOCK MODELS

# RELATIVE PARTITION RATES

clockRateID.a <- function(a) {
  paste0("clockRate.", get.gtree.PartitionRateM(a)$id)
}  
  
clockRateID.u <- function(u) {
  clocks <- get.gtree.PartitionRateMs()
  stopifnot(u <= length(clocks))
  paste0("clockRate.", clocks[[u]]$id)
}


# BRANCH RATE HETEROGENEITY

# TODO relaxed lognorm etc




######################################################
# OPERATORS
# STACEY

# STANDARD  for smcTree
nodeReheightID <-            function() { "Operator.nodeReheight" }
threeBranchAdjusterID <-     function() { "Operator.ThreeBranchAdjuster" }
nodesNudgeID <-              function() { "Operator.nodesNudge" }
coordinatedPruneRegraftID <- function() { "Operator.coordinatedPruneRegraft" }
focusedScalerID <-           function() { "Operator.focusedScaler" }
#heightsWarperID <-           function() { "Operator.heightsWarper" }

popSFScalerID <-             function(scale) { paste0("Operator.popSFScaler.",            scale) }
bdcGrowthRateScalerID <-     function(scale) { paste0("Operator.bdcGrowthScaler",         scale) }
bdcRelDeathScalerID <-       function(scale) { paste0("Operator.bdcReldeathScaler",       scale) }
bdcRelDeathRandWalkerID  <-  function(halfw) { paste0("Operator.bdcRelDeathRandWalker",   halfw) }
bdcCollapseWtScalerID <-     function(scale) { paste0("Operator.bdcCollapseWtScaler",     scale) }
bdcCollapseWtRandWalkerID <- function(halfw) { paste0("Operator.bdcCollapseWtRandWalker", halfw) }
bdcOriginHeightScalerID <-   function(scale) { paste0("Operator.bdcOriginHtScaler",       scale) }


# STANDARD  for gene trees
clockRateScalerID <- function(g) { paste0("Operator.clockRateScaler.", g) }
treeScalerID <- function(g) { paste0("Operator.treeScaler.", g) }
treeRootScalerID <- function(g) { paste0("Operator.treeRootScaler.", g) }
upClock.downHeightsID <- function(g) { paste0("Operator.upClock.downHeights.", g) }

subtreeSlideID <- function(g) { paste0("Operator.subtreeSlide.", g) }
treeScalerID <- function(g) { paste0("Operator.treeScaler.", g) }
wideID <- function(g) { paste0("Operator.wide.", g) }
narrowID <- function(g) { paste0("Operator.narrow.", g) }
WilsonBaldingID <- function(g) { paste0("Operator.WilsonBalding.", g) }
uniformID <- function(g) { paste0("Operator.uniform.", g) }

kappaScalerID.u <- function(u) { paste0("Operator.kappaScaler.", u) }
frequenciesExchangerID <- function(g) { paste0("Operator.frequenciesExchanger.", g) }


# STANDARD  for smcTree and gene trees
upGrowthClocks.downPopsHeightsID  <- function(scale) { paste0("Operator.upGrowthClocks.downPopsHeights", scale) }



#############################################################################
# LOGGERS

mainloggerID <- function() {
  "Logger.trace"
}


smctreeloggerID <- function() {
  "Logger.smcTree"
}


gtreeloggerID <- function(g) {
  paste0("Logger.gtree.", g)
}


screenloggerID <- function() {
  "Logger.screen"
}


