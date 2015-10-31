



#########################################################################
# for XML lib
#########################################################################


# general 

get.child <- function(node, name) {
	stopifnot(!is.null(node$children))
	stopifnot(is.numeric(node$children[name]))
	TheAnalysisStructureNodes[[ node$children[name] ]]
}

get.childvalue <- function(node, name) {
  as.vector(get.child(node, name)$value)
}


# analysis elements other than partitions
get.runoption <- function(name) {
	rops <- get.child(TheAnalysisStructureNodes[[1]], "run.options")
	get.childvalue(rops, name)
}



get.taxa.table <- function() {
	ttb.i <- index.from.kind.id("TaxaTable", TheAnalysisStructureNodes[[1]]$id)
	TheAnalysisStructureNodes[[ttb.i]]
}



get.minclusters <- function() {
	tt <- get.taxa.table()$value
	mcs <- tt[,2]
	unique(mcs)
}


# partitions - data


get.num.alignments <- function() {
	nrow(ThePartitionIndices)
}

get.data.fpaths <- function() {
	dp.i <- TheAnalysisStructureNodes[[1]]$children["data.dpath"]
	dpath <- TheAnalysisStructureNodes[[dp.i]]$value
	nFiles <- nrow(ThePartitionIndices)
	fpaths <- rep("", nFiles)
	for (a in 1:nFiles) {
		dfn.i <- ThePartitionIndices[a, "DataFileNames"]
		fname <- TheAnalysisStructureNodes[[dfn.i]]$value
		fpaths[a] <- paste0(dpath, "/", fname)
	}
	fpaths
}



get.gtreeprior.for.first.partition <- function() {
	gtp1.i <- ThePartitionIndices[1, "GTreePriors"]
	prr.i <- TheAnalysisStructureNodes[[gtp1.i]]$children["prior"]
	TheAnalysisStructureNodes[[prr.i]]
}


# gets branching (as in branching process) model for prior for first gene tree, 
# which, for STACEY, is the birth-death-collapse model for the SMC-tree. 
get.bdc.model <- function() {
	prr <- get.gtreeprior.for.first.partition()
	get.child(prr, "branching")
}




get.smc.coalescent <- function() {
	prr <- get.gtreeprior.for.first.partition()
	get.child(prr, "smccoal")
}



# These 11 functions return nodes of 11 kinds. 'a' is an alignment index. 

get.gtree.Partition <- function(a) {
	TheAnalysisStructureNodes[[ ThePartitionIndices[a, "Partitions"] ]]
}

get.gtree.DataFileName <- function(a) {
	TheAnalysisStructureNodes[[ ThePartitionIndices[a, "DataFileNames"] ]]
}

get.gtree.GTree <- function(a) {
	TheAnalysisStructureNodes[[ ThePartitionIndices[a, "GTrees"] ]]
}

get.gtree.Clock <- function(a) {
	TheAnalysisStructureNodes[[ ThePartitionIndices[a, "Clocks"] ]]
}

get.gtree.SiteM <- function(a) {
	TheAnalysisStructureNodes[[ ThePartitionIndices[a, "SiteMs"] ]]
}

get.gtree.GTreePrior <- function(a) {
	TheAnalysisStructureNodes[[ ThePartitionIndices[a, "GTreePriors"] ]]
}

get.gtree.BranchRM <- function(a) {
	TheAnalysisStructureNodes[[ ThePartitionIndices[a, "BranchRMs"] ]]
}

get.gtree.Ploidy <- function(a) {
	TheAnalysisStructureNodes[[ ThePartitionIndices[a, "Ploidys"] ]]
}

get.gtree.PartitionRateM <- function(a) {
	TheAnalysisStructureNodes[[ ThePartitionIndices[a, "PartitionRateMs"] ]]
}

get.gtree.SubstM <- function(a) {
	TheAnalysisStructureNodes[[ ThePartitionIndices[a, "SubstMs"] ]]
}

get.gtree.SiteHet <- function(a) {
	TheAnalysisStructureNodes[[ ThePartitionIndices[a, "SiteHets"] ]]
}


# These 11 functions return lists of nodes of 
# 11 different kinds. They are ones which exist
# per alignment, but are not necessarily unique
# to the alignment. The returned lists are lists of the
# unique ones which occur. They are parallel to
# get.Partitions() etc in analysis-structure-lib. Those are 
# for the user. These are for making XML. The difference is that these
# include which alignment (partition) is the first to use ir.

get.GTrees.with.first.partition <- function() {
	get.standard.nodes.with.first.partition("GTrees")
}

get.GTreePriors.with.first.partition <- function() {
	get.standard.nodes.with.first.partition("GTreePriors")
}

get.gtree.BranchRMs.with.first.partition <- function() {
	get.standard.nodes.with.first.partition("BranchRMs")
}


get.gtree.Ploidys.with.first.partition <- function() {
	get.standard.nodes.with.first.partition("Ploidys")
}


get.gtree.PartitionRateMs.with.first.partition <- function() {
	get.standard.nodes.with.first.partition("PartitionRateMs")
}


get.gtree.SiteMs.with.first.partition <- function() {
	get.standard.nodes.with.first.partition("SiteMs")
}
get.gtree.SubstMs.with.first.partition <- function() {
	get.standard.nodes.with.first.partition("SubstMs")
}

get.gtree.SiteHets.with.first.partition <- function() {
	get.standard.nodes.with.first.partition("SiteHets")
}



#############################  helper functions #############################


get.standard.nodes.with.first.partition <- function(kind) {
	nodes <- NULL
	nodeIDs <- character(0)
	for (a in 1:get.nof.alignments()) {
		node <- TheAnalysisStructureNodes[[ ThePartitionIndices[a, kind] ]]
		if (!(node$id %in% nodeIDs)) {
			nodeIDs <- c(nodeIDs, node$id)
			node.aug <- c(node, list(first.partition=a))
			nodes <- c(nodes, list(node.aug))
		}
	}
	nodes
}

