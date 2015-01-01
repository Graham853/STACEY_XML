#  TODO change a lot to use user-sipplied IDs

##################################################
# PARAMETERS
# STACEY 

smcTreeID <- function() { "smcTree" }
taxonSetOfSetsID <- function() { "taxonSetOfSets" }
# spp/min-cluster and sequences IDs not done here

bdcModelID <- function() { "BirthDeathCollapseModel" }
bdcGrowthID <- function() { "bdcGrowthRate.smcTree" }
bdcRelDeathID <- function() { "bdcRelativeDeathRate.smcTree" }
bdcCollapseWtID <- function() { "bdcCollapseWeight.smcTree" }
bdcOriginHtID <- function() { "bdcOriginHeight.smcTree" }
#bdcCollapseHtID <- function() { "collapseHeight" } an attribute, not element or parameter

popSFID <- function() { "popPriorScale" }

smcCoalescentID <- function() { "smcCoalescent" }
InvGammaComponentID <- function(c) { paste0("InvGammaComponent.", c) }



######################################################

alignmentID <- function(a) { paste0("algmt.", a) } 

# GENE TREES (may be linked across alignments)

geneTreeID <- function(g) { paste0("geneTree.", g) }
geneTreeCoalFactorID <- function(g) { paste0("gTreeCF.", g) }
geneTaxonSetID <- function(g) { paste0("geneTaxonSet.", g) }

geneTreeLhoodID <- function(g) { paste0("genetreeLhood.", g) }

# SITE HETEROGENEITY MODELS


# SUBSTITUTION MODELS

kappaID <- function(b) { 
  substs <- get.substs()
  stopifnot(b <= length(substs))
  paste0("kappa.", substs[[b]]$id)
}

frequenciesID <- function(b) {
  substs <- get.substs()
  stopifnot(b <= length(substs))
  paste0("frequencies.", substs[[b]]$id)
}

# CLOCK MODELS

clockRateID <- function(c) {
  clocks <- get.clocks()
  stopifnot(c <= length(clocks))
  paste0("clockRate.", clocks[[c]])
}


######################################################
# OPERATORS
# STACEY

nodesNudgeID <- function() { "nodesNudge" }
coordinatedPruneRegraftID <- function() { "coordinatedPruneRegraft" }
focusedScalerID <- function() { "focusedScaler" }
popSFScalerID <- function() { "popSFScaler" }
bdcGrowthScalerID <- function() { "bdcGrowthScaler" }
bdcReldeathScalerID <- function() { "bdcReldeathScaler" }
bdcCollapseWtScalerID <- function() { "bdcCollapseWtScaler" }
bdcOriginHtScalerID <- function() { "bdcOriginHtScaler" }


# STANDARD: for smcTree
nodeReheightID <- function() { "nodeReheight" }

# STANDARD: for gene trees
clockRateScalerID <- function(g) { paste0("clockRateScaler.", g) }
treeScalerID <- function(g) { paste0("treeScaler.", g) }
treeRootScalerID <- function(g) { paste0("treeRootScaler.", g) }
upClock.downHeightsID <- function(g) { paste0("upClock.downHeights.", g) }

subtreeSlideID <- function(g) { paste0("subtreeSlide.", g) }
treeScalerID <- function(g) { paste0("treeScaler.", g) }
wideID <- function(g) { paste0("wide.", g) }
narrowID <- function(g) { paste0("narrow.", g) }
WilsonBaldingID <- function(g) { paste0("WilsonBalding.", g) }
uniformID <- function(g) { paste0("uniform.", g) }

kappaScalerID <- function(g) { paste0("kappaScaler.", g) }
frequenciesExchangerID <- function(g) { paste0("frequenciesExchanger.", g) }


# STANDARD: for smcTree and gene trees
upGrowthClocks.downPopsHeightsID  <- function() { "upGrowthClocks.downPopsHeights" }


IDtoREF <- function(ID) {
  paste0("@", ID)
}








