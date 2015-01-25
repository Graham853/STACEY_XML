
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
  fname <- TheAnalysisStructure$alignment.table$alignments[[a]]$file
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
  paste0("Partition:", id)
}

# data for partition a
partitiondataID.a <- function(a) {
  paste0("Data:", dataID.a(a))
}

# data for first partition belonging gene tree g
partitiondataID.g <- function(g) {
  a <- get.gtrees()[[g]]$first.partition
  paste0("Data:", dataID.a(a))
}



sequenceID.a <- function(a, taxonname) {
  paste0("Seq:", dataID.a(a), ".", taxonname)
}





initialsmctreeID <- function() {
  "init:smcTree"
}


initialsmctreepopmodelID <- function() {
  "init:smcTree.popmodel"
}


initialsmctreepopsizeID <- function() {
  "init:smcTree:popsize"
}


initialgtreeID <- function(g) {
  paste0("init:gTree", g)
}


initialgtreepopmodelID <- function(g) {
  paste0("init:gTree.popmodel", g)
}


initialgtreepopsizeID <- function(g) {
  paste0("init:gTree.popsize", g)
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
  smct <- prior.of.first.gtree()
  xmlIDFromElement(smct)
}

taxonSetOfSetsID <- function() { "taxonSetOfSets" }
# spp/min-cluster and sequences IDs not done here

bdcModelID <- function() {
  smct <- prior.of.first.gtree()
  xmlIDFromElement(smct$branching)
}


bdcGrowthID <- function() { 
  smct <- prior.of.first.gtree()
  xmlIDfromKindID(smct$branching$kind, smct$branching$growthrate$id)
}


bdcRelDeathID <- function() {
  smct <- prior.of.first.gtree()
  xmlIDfromKindID(smct$branching$kind, smct$branching$reldeath$id)
}


bdcCollapseWtID <- function() { 
  smct <- prior.of.first.gtree()
  xmlIDfromKindID(smct$branching$kind, smct$branching$w$id)
}


bdcOriginHtID <- function() {
  smct <- prior.of.first.gtree()
  xmlIDfromKindID(smct$branching$kind, smct$branching$ohID)
}



smcCoalescentID <- function() {
  smct <- prior.of.first.gtree()
  xmlIDFromElement(smct$smccoal) 
}


InvGammaComponentID <- function(c) { 
  smct <- prior.of.first.gtree()
  paste0(smct$smccoal$kind, ".", smct$smccoal$invgammamix$id, ".", c) 
}


popSFID <- function() {
  smct <- prior.of.first.gtree()
  xmlIDfromKindID(smct$smccoal$kind, smct$smccoal$popSF$id) 
}



######################################################


########################################################

# GENE TREES. 




geneTreeID.a <- function(a) {
  gtree <- TheAnalysisStructure$alignment.table$alignments[[a]]$gtree
  xmlIDFromElement(gtree)
}

geneTreeID.g <- function(g) {
  gtree <- get.gtrees()[[g]]
  xmlIDFromElement(gtree)
}



geneTreeBranchRM <- function(g) {
  gtree <- get.gtrees()[[g]]
  xmlIDFromElement(gtree$branchRM$model)
}



geneTreeLhoodID <- function(g) {
  gtree <- get.gtrees()[[g]]
  paste0("GTreeLhood:", gtree$id)
}


geneTreeCoalFactorID.g <- function(g) {
  gtree <- get.gtrees()[[g]]
  paste0("GTreePloidy:", gtree$id)
}


# referred to in init to make random trees
geneTaxonSetID.g <- function(g) {
  gtree <- get.gtrees()[[g]]
  paste0("GTreeTaxonSet:", gtree$id)
}



####################################################
##### SITE MODELS


sitemodelID.a  <- function(a) {
  siteM <- TheAnalysisStructure$alignment.table$alignments[[a]]$siteM
  xmlIDFromElement(siteM)
}

sitemodelID.u <- function(u) {
  siteM <- get.siteMs()[[u]]
  xmlIDFromElement(siteM)
}

# SITE HET

sitehetMuID.a <- function(a) {
  sitehet <- TheAnalysisStructure$alignment.table$alignments[[a]]$siteM$sitehet
  stopifnot(is.character(sitehet$model))
  if (sitehet$model == "None") {
    paste0(sitehet$kind, ".mu.", sitehet$id)
  } else {
    #TODO 
  }  
}


sitehetGammaShapeID.a <- function(a) {
  sitehet <- TheAnalysisStructure$alignment.table$alignments[[a]]$siteM$sitehet
  stopifnot(is.character(sitehet$model))
  if (sitehet$model == "None") {
    paste0(sitehet$kind, ".GammaShape.", sitehet$id)
  } else {
    #TODO 
  }  
}

sitehetPropInvID.a <- function(a) {
  sitehet <- TheAnalysisStructure$alignment.table$alignments[[a]]$siteM$sitehet
  stopifnot(is.character(sitehet$model))
  if (sitehet$model == "None") {
    paste0(sitehet$kind, ".PropInv.", sitehet$id)
  } else {
    #TODO 
  }  
}




# SUBSTITUTION MODELS 

substsmodelID.a <- function(a) {
  model <- TheAnalysisStructure$alignment.table$alignments[[a]]$siteM$subst$model
  xmlIDFromElement(model)
}

substsmodelID.u <- function(u) {
  model <- get.siteMs()[[u]]$subst$model
  xmlIDFromElement(model)
}


kappaID.a <- function(a) {
  model <- TheAnalysisStructure$alignment.table$alignments[[a]]$siteM$subst$model
  xmlIDfromKindID(model$kind, model$kappa$id)
}

kappaID.u <- function(u) { 
  model <- get.siteMs()[[u]]$subst$model
  xmlIDfromKindID(model$kind, model$kappa$id)
}


# there are two things that are 'frequencies'. model, param, i think
# <frequencies frequencies="@freqParameter.s:U2SQg5t10r1-1" id="estimatedFreqs.s:U2SQg5t10r1-1" spec="Frequencies"/>


frequenciesModelID.a <- function(a) {
  model <- TheAnalysisStructure$alignment.table$alignments[[a]]$siteM$subst$model
  paste0("FrequenciesModel:", model$freqs$id)
}


frequenciesParamID.a <- function(a) {
  model <- TheAnalysisStructure$alignment.table$alignments[[a]]$siteM$subst$model
  xmlIDfromKindID(model$kind, model$freqs$id)
}

frequenciesParamID.u <- function(u) {
  model <- get.siteMs()[[u]]$subst$model
  xmlIDfromKindID(model$kind, model$freqs$id)
}



# CLOCK RATES (RELATIVE PARTITION RATES)

clockRateID.a <- function(a) {
  paste0("clockRate:", TheAnalysisStructure$alignment.table$alignments[[a]]$clock$id)
}  
  
clockRateID.c <- function(c) {
  clocks <- get.clocks()
  stopifnot(c <= length(clocks))
  paste0("clockRate:", clocks[[c]]$id)
}




######################################################
# OPERATORS
# STACEY

nodesNudgeID <- function() { "Operator:nodesNudge" }
coordinatedPruneRegraftID <- function() { "Operator:coordinatedPruneRegraft" }
focusedScalerID <- function() { "Operator:focusedScaler" }
popSFScalerID <- function() { "Operator:popSFScaler" }
bdcGrowthScalerID <- function() { "Operator:bdcGrowthScaler" }
bdcReldeathScalerID <- function() { "Operator:bdcReldeathScaler" }
bdcCollapseWtScalerID <- function() { "Operator:bdcCollapseWtScaler" }
bdcOriginHtScalerID <- function() { "Operator:bdcOriginHtScaler" }


# STANDARD: for smcTree
nodeReheightID <- function() { "Operator:nodeReheight" }

# STANDARD: for gene trees
clockRateScalerID <- function(g) { paste0("Operator:clockRateScaler.", g) }
treeScalerID <- function(g) { paste0("Operator:treeScaler.", g) }
treeRootScalerID <- function(g) { paste0("Operator:treeRootScaler.", g) }
upClock.downHeightsID <- function(g) { paste0("Operator:upClock.downHeights.", g) }

subtreeSlideID <- function(g) { paste0("Operator:subtreeSlide.", g) }
treeScalerID <- function(g) { paste0("Operator:treeScaler.", g) }
wideID <- function(g) { paste0("Operator:wide.", g) }
narrowID <- function(g) { paste0("Operator:narrow.", g) }
WilsonBaldingID <- function(g) { paste0("Operator:WilsonBalding.", g) }
uniformID <- function(g) { paste0("Operator:uniform.", g) }

kappaScalerID.u <- function(u) { paste0("Operator:kappaScaler.", u) }
frequenciesExchangerID <- function(g) { paste0("Operator:frequenciesExchanger.", g) }


# STANDARD: for smcTree and gene trees
upGrowthClocks.downPopsHeightsID  <- function() { "Operator:upGrowthClocks.downPopsHeights" }



#############################################################################


mainloggerID <- function() {
  "Logger:trace"
}


smctreeloggerID <- function() {
  "Logger:smcTree"
}


gtreeloggerID <- function(g) {
  paste0("Logger:gtree.", g)
}


screenloggerID <- function() {
  "Logger:screen"
}


