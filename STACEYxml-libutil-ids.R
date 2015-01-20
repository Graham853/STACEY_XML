#  TODO change a lot to use user-sipplied IDs

##################################################
# PARAMETERS
# STACEY 

prior.of.first.gtree <- function() {
  TheAnalysisStructure$alignment.table$alignments[[1]]$gtree$prior
}


smcTreeID <- function() {
  smct <- prior.of.first.gtree()
  paste0(smct$kind, ".", smct$id)
}

taxonSetOfSetsID <- function() { "taxonSetOfSets" }
# spp/min-cluster and sequences IDs not done here

bdcModelID <- function() {
  smct <- prior.of.first.gtree()
  paste0(smct$branching$kind, ".", smct$branching$id)
}


bdcGrowthID <- function() { 
  smct <- prior.of.first.gtree()
  paste0(smct$branching$kind, ".", smct$branching$growthrate$id)
}


bdcRelDeathID <- function() {
  smct <- prior.of.first.gtree()
  paste0(smct$branching$kind, ".", smct$branching$reldeath$id)
}


bdcCollapseWtID <- function() { 
  smct <- prior.of.first.gtree()
  paste0(smct$branching$kind, ".", smct$branching$w$id)
}


bdcOriginHtID <- function() {
  smct <- prior.of.first.gtree()
  paste0(smct$branching$kind, ".", smct$branching$ohID)
}



smcCoalescentID <- function() {
  smct <- prior.of.first.gtree()
  paste0(smct$smccoal$kind, ".", smct$smccoal$id) 
}


InvGammaComponentID <- function(c) { 
  smct <- prior.of.first.gtree()
  paste0(smct$smccoal$kind, ".", smct$smccoal$invgammamix$id, ".", c) 
}


popSFID <- function() {
  smct <- prior.of.first.gtree()
  paste0(smct$smccoal$kind, ".", smct$smccoal$popSF$id) 
}



######################################################


########################################################

# GENE TREES. 

# data for first partition belonging to this gene tree
gtreedataID.g <- function(g) {
  a <- get.gtrees()[[g]]$first.partition
  paste0("Data", ".", TheAnalysisStructure$alignment.table$alignments[[a]]$file)
}


geneTreeID.a <- function(a) {
  gtree <- TheAnalysisStructure$alignment.table$alignments[[a]]$gtree
  paste0(gtree$kind, gtree$id)
}

geneTreeID.g <- function(g) {
  gtree <- get.gtrees()[[g]]
  paste0(gtree$kind, ".", gtree$id)
}


# 
# branchRM:List of 3
# ..$ kind : chr "BranchRM"
# ..$ id   : chr "1"
# ..$ model:List of 3
# .. ..$ kind: chr "StrictClock"
# .. ..$ id  : chr "1"
# .. ..$ rate: num 1


geneTreeLhoodID <- function(g) {
  gtree <- get.gtrees()[[g]]
  paste0("Lhood.", gtree$id)
}





geneTreeCoalFactorID.g <- function(g) {
  gtree <- get.gtrees()[[g]]
  paste0("gTreeCF.", gtree$id)
}



geneTaxonSetID.g <- function(g) {
  gtree <- get.gtrees()[[g]]
  paste0("geneTaxonSet.", gtree$id)
}





# TODO SITE HET


# SUBSTITUTION MODELS substsID.u(u) is u'th ID in alignment table

substsID.u <- function(u) {
  substs <- get.siteMs()
  stopifnot(u <= length(substs))
  paste0(substs[[u]]$kind, ".", substs[[u]]$id)
}


kappaID.a <- function(a) {
  subst <- TheAnalysisStructure$alignment.table$alignments[[a]]$siteM$subst
  paste0(subst$model$kind, subst$model$kappa$id)
}

kappaID.u <- function(u) { 
  subst <- get.siteMs()[[u]]$subst
  paste0(subst$model$kind, subst$model$kappa$id)
}

frequenciesID.a <- function(a) {
  subst <- TheAnalysisStructure$alignment.table$alignments[[a]]$siteM$subst
  paste0(subst$model$kind, subst$model$freqs$id)
}

frequenciesID.u <- function(u) {
  subst <- get.siteMs()[[u]]
  paste0(subst$model$kind, subst$model$freqs$id)
}



# CLOCK MODELS clocksID.c(c) is c'th ID in alignment table

clocksID.c <- function(c) {
  clocks <- get.clocks()
  stopifnot(c <= length(clocks))
  clocks[[c]]$id
}
  
clockRateID.a <- function(a) {
  paste0("frequencies.", TheAnalysisStructure$alignment.table$alignments[[a]]$clock$id)
}  
  
clockRateID.c <- function(c) {
  clocks <- get.clocks()
  stopifnot(c <= length(clocks))
  paste0("clockRate.", clocks[[c]]$id)
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

kappaScalerID.u <- function(u) { paste0("kappaScaler.", u) }
frequenciesExchangerID <- function(g) { paste0("frequenciesExchanger.", g) }


# STANDARD: for smcTree and gene trees
upGrowthClocks.downPopsHeightsID  <- function() { "upGrowthClocks.downPopsHeights" }


IDtoREF <- function(ID) {
  paste0("@", ID)
}








