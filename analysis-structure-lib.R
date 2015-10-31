# This provides lots of functions to build an analysis structure.
library(hash)
library(stringr)
library(ape)


# July 2015. Rewrite to use a list of nodes, not a list of lists of lists...

# A node is a list with elements:   kind, id, link, value, children

# kind and id are strings. (type character, length 1.) The other three take
# values in three ways
#             link     value     children
#             index    NULL      NULL
#             NULL     something NULL
#             NULL     NULL      indicies
#             NULL     NULL      NULL
# 'index' and 'indices' are integers specifying other nodes. 
# 'something' is a vector (often length 1)of type "logical", "integer", "numeric", "complex", 
# "character" or "raw". (is.atomic(), but not NULL. I don't imagine raw being used.)
# link, value, children are always present and explicitly NULL as appropriate, not omitted.
# All NULL means an incomplete internal node.


TheAnalysisStructureNodes <- NULL
ThePartitionIndices <- NULL
TheNodeHashTable <- hash()

TheAnalysis <- function(id, data.dpath, alignment.table, taxa.table, run.options) {
  stopifnot(is.character(data.dpath))
  stopifnot(length(data.dpath)==1)
  stopifnot(is.character(data.fnames))
  stopifnot(is.character(run.options))
  taxa.table <- as.matrix(taxa.table)
  stopifnot(ncol(taxa.table)==2)
  
  # the root, with 4 children. Data directory and taxa table are leaves.
  tas.i <- add.new.node("TheAnalysisStructure", id)
  ddr.i <- add.leaf("DataDirectory", id, data.dpath)
  atb.i <- add.new.node("AlignmentTable", id)
  ttb.i <- add.leaf("TaxaTable", id, taxa.table)
  rop.i <- add.new.node("RunOptions", id)
  children <- c(ddr.i, atb.i, ttb.i, rop.i)
  names(children) <- c("data.dpath", "alignment.table", "taxa.table", "run.options")
  attach.children(tas.i, children)
  
  # The run options. 9 strings
  stopifnot(is.character(run.options["sampledgtrees.fpathbase"]))
  stopifnot(is.character(run.options["sampledsmctrees.fpath"]))
  stopifnot(is.character(run.options["sampledparams.fpath"]))
  stopifnot(is.character(run.options["chainlength"]))
  stopifnot(is.character(run.options["store.every"]))
  stopifnot(is.character(run.options["params.logevery"]))
  stopifnot(is.character(run.options["smctree.logevery"]))
  stopifnot(is.character(run.options["gtrees.logevery"]))
  stopifnot(is.character(run.options["screen.logevery"]))
  gts.fpb <- add.leaf("sampledgtrees.fpathbase", id, run.options["sampledgtrees.fpathbase"])
  sts.fp  <- add.leaf("sampledsmctrees.fpath", id, run.options["sampledsmctrees.fpath"])
  spm.fp  <- add.leaf("sampledparams.fpath", id, run.options["sampledparams.fpath"])
  ch.len  <- add.leaf("chainlength", id, run.options["chainlength"])
  store.evry <- add.leaf("store.every", id, run.options["store.every"])
  pm.evry <- add.leaf("params.logevery", id, run.options["params.logevery"])
  smct.evry <- add.leaf("smctree.logevery", id, run.options["smctree.logevery"])
  gt.evry <- add.leaf("gtrees.logevery", id, run.options["gtrees.logevery"])
  scrn.evry <- add.leaf("screen.logevery", id, run.options["screen.logevery"])
  children <- c(gts.fpb, sts.fp, spm.fp, ch.len, store.evry, pm.evry, smct.evry, gt.evry, scrn.evry)
  names(children) <- c("sampledgtrees.fpathbase", "sampledsmctrees.fpath", "sampledparams.fpath",
  										 "chainlength", "store.every",
  										 "params.logevery", "smctree.logevery", "gtrees.logevery", "screen.logevery")
  attach.children(rop.i, children)
  
  # The alignment, the main part.
  # For each data file, add a partition with 4 children. File names are leaves
  
  ThePartitionIndices <<- matrix(-1, ncol=11, nrow=length(data.fnames))
  colnames(ThePartitionIndices) <<- c("Partitions", "DataFileNames", "GTrees", "Clocks", "SiteMs", 
                                      "GTreePriors", "BranchRMs", "Ploidys",
                                      "PartitionRateMs", 
                                      "SubstMs", "SiteHets")
  
  for (a in 1:length(data.fnames)) {
    id <- names(data.fnames)[a]
    ThePartitionIndices[a, "Partitions"]    <<- ptn.i <- add.new.node("Partition", id)
    ThePartitionIndices[a, "DataFileNames"] <<- dfn.i <- add.leaf("DataFileName", id, data.fnames[a])
    ThePartitionIndices[a, "GTrees"]        <<- gtr.i <- add.new.node("GTree", id)
    ThePartitionIndices[a, "Clocks"]        <<- clk.i <- add.new.node("Clock", id)
    ThePartitionIndices[a, "SiteMs"]        <<- stm.i <- add.new.node("SiteM", id)
    children <- c(dfn.i, gtr.i, clk.i, stm.i)
    names(children) <- c("DataFileNames", "GTrees", "Clocks", "SiteMs")
    attach.children(ptn.i, children)
  }
  children <- ThePartitionIndices[ , "Partitions"]
  names <- rep("", length(children)) 
  for (i in 1:length(children)) {
  	names[i] <- TheAnalysisStructureNodes[[children[i]]]$id
  }
  names(children) <- names
  attach.children(atb.i, children)
 
  
  # For each data file, the gene tree has three children to be filled later
  for (a in 1:length(data.fnames)) {
    gtr.i <- ThePartitionIndices[a, "GTrees"]
    ThePartitionIndices[a, "GTreePriors"] <<- gpr.i <- add.new.node("GTreePrior", names(data.fnames)[a])
    ThePartitionIndices[a, "BranchRMs"]   <<- gbr.i <- add.new.node("BranchRM", names(data.fnames)[a])
    ThePartitionIndices[a, "Ploidys"]     <<- gpy.i <- add.new.node("Ploidy", names(data.fnames)[a])
    children <- c(gpr.i, gbr.i, gpy.i)
    names(children) <- c("GTreePrior", "BranchRM", "Ploidy")
    attach.children(gtr.i, children)
  }

  
  # For each data file, the clock model has a child to be filled later
  for (a in 1:length(data.fnames)) {
    clk.i <- ThePartitionIndices[a, "Clocks"]
    ThePartitionIndices[a, "PartitionRateMs"] <<- prm.i <- add.new.node("PartitionRateM", names(data.fnames)[a])
    children <- prm.i
    names(children) <- "PartitionRateM"
    attach.children(clk.i, children)
  }

  # For each data file, the site model has two children to be filled later
  for (a in 1:length(data.fnames)) {
    stm.i <- ThePartitionIndices[a, "SiteMs"]
    ThePartitionIndices[a, "SubstMs"]  <<- gsu.i <- add.new.node("SubstM", names(data.fnames)[a])
    ThePartitionIndices[a, "SiteHets"] <<- ghe.i <- add.new.node("SiteHet", names(data.fnames)[a])
    children <- c(gsu.i, ghe.i)
    names(children) <- c("SubstM", "SiteHet") 
    attach.children(stm.i, children)
  }
}



##################################################################


get.nof.alignments <- function() {
	nrow(ThePartitionIndices)
}

# these 11 functions return lists of nodes of 
# 11 different kinds. They are ones which exist
# per alignment, but re not necessarily unique
# to the alignment. The returned lists are lists of the
# unique ones which occur.

get.Partitions <- function() {
	get.standard.nodes("Partitions")
}

get.DataFileNames <- function() {
	get.standard.nodes("DataFileNames")
}

get.GTrees <- function() {
	get.standard.nodes("GTrees")
}

get.Clocks <- function() {
	get.standard.nodes("Clocks")
}

get.SiteMs <- function() {
	get.standard.nodes("SiteMs")
}

get.GTreePriors <- function() {
	get.standard.nodes("GTreePriors")
}

get.gtree.BranchRMs <- function() {
	get.standard.nodes("BranchRMs")
}

get.gtree.Ploidys <- function() {
	get.standard.nodes("Ploidys")
}

get.gtree.PartitionRateMs <- function() {
	get.standard.nodes("PartitionRateMs")
}

get.gtree.SubstMs <- function() {
	get.standard.nodes("SubstMs")
}

get.gtree.SiteHets <- function() {
	get.standard.nodes("SiteHets")
}


######################################################################
################### Building functions ###############################
######################################################################


##########################################################################
# Gene trees and their prior(s)

GTreePrior <- function(id, prior) { 
	x <- list(kind="GTreePrior", id=id, prior=prior)
	if (bbas.real.kids(x)) {
		gtp.i <- index.from.kind.id("GTreePrior", id)
		prr.i <- add.new.node(prior$kind, prior$id)
		children <- prr.i
		names(children) <- names(x[c(-1,-2)])
		attach.children(gtp.i, children)		
	}
}


SMCtree <- function(id, branching=NULL, smccoal=NULL) {
  x <- list(kind="SMCtree", id=id, branching=branching, smccoal=smccoal)
  if (bbas.real.kids(x)) {
    sct.i <- index.from.kind.id("SMCtree", id)
    bra.i <- add.new.node(branching$kind, branching$id)
    scl.i <- add.new.node(smccoal$kind, smccoal$id)
    children <- c(bra.i, scl.i)
    names(children) <- names(x[c(-1,-2)])
    attach.children(sct.i, children)
  }
  invisible(x)
}



BDCPrior <- function(id, growthrate=NULL, reldeath=NULL, w=NULL, oh=NULL, eps=NULL) {
  x <- list(kind="BDCPrior", id=id, growthrate=growthrate, reldeath=reldeath, w=w, oh=oh, eps=eps)
  if (bbas.real.kids(x)) {
    bdc.i <- index.from.kind.id("BDCPrior", id)
    grw.i <- add.new.node(growthrate$kind, growthrate$id)
    rdh.i <- add.new.node(reldeath$kind, reldeath$id)
    cwt.i <- add.new.node(w$kind, w$id)
    oht.i <- add.leaf("SMCTreeOriginHeight", oh, 0.05)
    eps.i <- add.leaf("BDCollapseHeight", "eps", eps) 
    children <- c(grw.i, rdh.i, cwt.i, oht.i, eps.i)
    names(children) <- names(x[c(-1,-2)])
    attach.children(bdc.i, children)
  }
  invisible(x)
}  


SMCCoalescent <- function(id, invgammamix=NULL, popSF=NULL) {
  x <- list(kind="SMCCoalescent", id=id, invgammamix=invgammamix, popSF=popSF)
  if (bbas.real.kids(x)) {
    sct.i <- index.from.kind.id("SMCCoalescent", id)
    ivg.i <- add.new.node(invgammamix$kind, invgammamix$id)
    psf.i <- add.new.node(popSF$kind, popSF$id)
    children <- c(ivg.i, psf.i)
    names(children) <- names(x[c(-1,-2)])
    attach.children(sct.i, children)
  }
  invisible(x)
}    


##########################################################################
# Clock models for branch rate heterogeneity

BranchRM <- function(id, brhetM=NULL) {
  x <- list(kind="BranchRM", id=id, brhetM=brhetM)
  if (bbas.real.kids(x)) {
    bra.i <- index.from.kind.id("BranchRM", id)
    bhe.i <- add.new.node(brhetM$kind, brhetM$id)
    children <- bhe.i
    names(children) <- names(x[c(-1,-2)])
    attach.children(bra.i, children)  
  }
  invisible(x)
}


# a model for branch rate heterogeneity
StrictClock <- function(id, rate=NULL) {
  x <- list(kind="StrictClock", id=id, rate=rate)
  if (bbas.real.kids(x)) {
    stopifnot(is.numeric(rate))
    sck.i <- index.from.kind.id("StrictClock", id)
    scr.i <- add.leaf("StrictClockRate", id, rate)
    children <- scr.i
    names(children) <- names(x[c(-1,-2)])    
    attach.children(sck.i, children) 
  }
  invisible(x)
}


##########################################################################
# Clock models for relative partition rates 


PartitionRateM <- function(id, rate=NULL) {
  x <- list(kind="PartitionRateM", id=id, rate=rate)
  if (bbas.real.kids(x)) {
    prm.i <- index.from.kind.id("PartitionRateM", id)
    if (is.numeric(rate)) {
      ptr.i <- add.leaf("Fixed", id, rate)
    } else {
      ptr.i <- add.new.node(rate$kind, rate$id)
    }
    children <- ptr.i
    names(children) <- names(x[c(-1,-2)])    
    attach.children(prm.i, children)  
  }
  invisible(x)
}


############################################################################
# ploidy

Ploidy  <- function(id, value=NULL) {
  x <- list(kind="Ploidy", id=id, value=value)
  if (bbas.real.kids(x)) {
    pdy.i <- index.from.kind.id("Ploidy", id)
    stopifnot(is.numeric(value))
    pdv.i <- add.leaf("PloidyValue", id, value)
    children <- pdv.i
    names(children) <- names(x[c(-1,-2)])    
    attach.children(pdy.i, children)  
  }
  invisible(x)
}



##########################################################################
# site models: substitution and site rate heterogeneity models


SubstM <- function(id, model=NULL) {
  x <- list(kind="SubstM", id=id, model=model)
  if (bbas.real.kids(x)) {
    sub.i <- index.from.kind.id("SubstM", id) 
    smd.i <- add.new.node(model$kind, model$id)
    children <- smd.i
    names(children) <- names(x[c(-1,-2)])    
    attach.children(sub.i, children)   
  }
  invisible(x)
}


HKY <- function(id, freqs=NULL, kappa=NULL) {
  x <- list(kind="HKY", id=id, freqs=freqs, kappa=kappa)
  if (bbas.real.kids(x)) {
    smd.i <- index.from.kind.id("HKY", id) 
    frq.i <- add.new.node(freqs$kind, freqs$id)
    kpa.i <- add.new.node(kappa$kind, kappa$id)
    children <- c(frq.i, kpa.i)
    names(children) <- names(x[c(-1,-2)])    
    attach.children(smd.i, children)   
  }
  invisible(x)
}


FixedFrequencies <- function(id, howfixed=NULL) {
  x <- list(kind="FixedFrequencies", id=id, howfixed=howfixed)
  if (bbas.real.kids(x)) {
    fill.in.value("FixedFrequencies", id, howfixed)
  }
  invisible(x)  
}


SiteHet <- function(id, model=NULL) {
  x <- list(kind="SiteHet", id=id, model=model)
  if (bbas.real.kids(x)) {
    if (is.character(model)) {
      fill.in.value("SiteHet", id, model)
    }
  }
  invisible(x)
}


##########################################################################
# priors, that is, distributions with all parameters specified.
# Mostly they are one-dimensional.
# These could have versions with hyper-params. 


Fixed <- function(id, value=NULL) {
  x <- list(kind="Fixed", id=id, value=value)
  if (bbas.real.kids(x)) {
    stopifnot(is.numeric(value)) 
    fill.in.value("Fixed", id, value)
  }
  invisible(x)
}


InvGammaMix <- function(id, weights=NULL, alphas=NULL, betas=NULL) {
  x <- list(kind="InvGammaMix", id=id, weights=weights, alphas=alphas, betas=betas)
  if (bbas.real.kids(x)) {
    n <- length(weights)
    if (n != length(alphas)  ||  n != length(betas)) {
      stop("Vectors weights, alphas, and betas must have same dimension in InvGammaMix.")
    }
    stopifnot(is.numeric(weights))
    stopifnot(is.numeric(alphas))
    stopifnot(is.numeric(betas))
    initval <- sum(weights*(betas/alphas)) / sum(weights)
    # mean is b/(a-1) for a>1. mode is b/(a+1). mm is between
    igm.i <- index.from.kind.id("InvGammaMix", id)
    wts.i <- add.leaf("InvGammaMixWeights", id, weights) 
    alp.i <- add.leaf("InvGammaMixAlphas", id, alphas) 
    bet.i <- add.leaf("InvGammaMixBetas", id, betas) 
    inv.i <- add.leaf("InvGammaMixInitValue", id, initval) 
    children <- c(wts.i, alp.i, bet.i, inv.i)
    names(children) <- c("weights", "alphas", "betas", "initval")
    attach.children(igm.i, children)
  }
  invisible(x)
}  




MeanOneGamma <- function(id, alpha=NULL) {
  x <- list(kind="MeanOneGamma", id=id, alpha=alpha)
  if (bbas.real.kids(x)) {
    initval <- 1.0 
    stopifnot(is.numeric(alpha)) # TODO hyper-priors
    mog.i <- index.from.kind.id("MeanOneGamma", id)
    alp.i <- add.leaf("MeanOneGammaAlpha", id, alpha) 
    inv.i <- add.leaf("MeanOneGammaInitValue", id, initval) 
    children <- c(alp.i, inv.i)
    names(children) <- c("alpha", "initval")
    attach.children(lgn.i, children) 
  }
  invisible(x)
} 




LogNorm <- function(id, meanlog=NULL, sdlog=NULL) {
  x <- list(kind="LogNorm", id=id, meanlog=meanlog, sdlog=sdlog)
  if (bbas.real.kids(x)) {
    stopifnot(is.numeric(meanlog))
    stopifnot(is.numeric(sdlog)) # TODO hyper-priors
    initval <- exp(meanlog)
    lgn.i <- index.from.kind.id("LogNorm", id)
    mlg.i <- add.leaf("LogNormMeanLog", id, meanlog) 
    slg.i <- add.leaf("LogNormSDLog", id, sdlog) 
    inv.i <- add.leaf("LogNormInitValue", id, initval) 
    children <- c(mlg.i, slg.i, inv.i)
    names(children) <- c("meanlog", "sdlog", "initval")
    attach.children(lgn.i, children) 
  }
  invisible(x)
}


Beta <- function(id, a=NULL, b=NULL) {
  x <- list(kind="Beta", id=id, a=a, b=b)
  if (bbas.real.kids(x)) {
    stopifnot(is.numeric(a))
    stopifnot(is.numeric(b)) # TODO hyper-priors
    initval <- a/(a+b)
    beta.i <- index.from.kind.id("Beta", id)
    a.i <- add.leaf("BetaA", id, a) 
    b.i <- add.leaf("BetaB", id, b) 
    inv.i <- add.leaf("BetaInitValue", id, initval)  
    children <- c(a.i, b.i, inv.i)
    names(children) <- c("a", "b", "initval")
    attach.children(beta.i, children)       
  }
  invisible(x)
}


Uniform <- function(id, lower=NULL, upper=NULL) {
  x <- list(kind="Uniform", id=id, lower=lower, upper=upper)
  if (bbas.real.kids(x)) {
    stopifnot(is.numeric(lower))
    stopifnot(is.numeric(upper)) # TODO hyper-priors
    initval <- (lower+upper)/2
    unf.i <- index.from.kind.id("Uniform", id)
    low.i <- add.leaf("UniformLower", id, lower)
    upp.i <- add.leaf("UniformUpper", id, upper)
    inv.i <- add.leaf("UniformInitValue", id, initval) 
    children <- c(low.i, upp.i, inv.i)
    names(children) <- c("lower", "upper", "initval")
    attach.children(unf.i, children)
  }
  invisible(x)
}



UniformUnitSimplex <- function(id, N=NULL) {
  x <- list(kind="UniformUnitSimplex", id=id, N=N)
  if (bbas.real.kids(x)) {
    initval <- rep(1/N, N) # TODO I think BEAST only understands single values, so rep() is redundant
    unf.i <- index.from.kind.id("UniformUnitSimplex", id)
    inv.i <- add.leaf("UniformUnitSimplexInitValue", id, initval)
    children <- inv.i
    names(children) <- "initval"
    attach.children(unf.i, children) 
  }
  invisible(x)
}



######################################################################
######################## utility functions   #########################
######################################################################



tidyup <- function() {
  clear(TheNodeHashTable)
}


cat.structure <- function(file=NULL) {
  cat("Analysis Structure as tree with links\n", file=file)
  cat.substructure(1, 0, file)
}


cat.substructure <- function(n, level, file) {
  if (level > 0) {
    for (i in 1:level) { cat("    ", sep="", file=file, append=TRUE) }
  }
  node <- TheAnalysisStructureNodes[[n]]
  cat.node(n, file)
  if (!is.null(node$children)) {
    for (c in 1:length(node$children)) {
      cat.substructure(node$children[c], level+1, file)
    }
  }
}



cat.node <- function(n, file) { 
  node <- TheAnalysisStructureNodes[[n]] 
  cat(sprintf("[[%d]]  %s  %s ", n, node$kind, node$id), file=file, append=TRUE)
  if (!is.null(node$link)) {
    cat("link", node$link, "\n", file=file, append=TRUE)
  } else if (!is.null(node$value)) {
    cat("value", node$value, "\n", file=file, append=TRUE)
  } else if (!is.null(node$children)) {
  	cat(" children:", file=file, append=TRUE)
  	for (c in 1:length(node$children)) {
  		cat(" ", names(node$children)[c], "=", node$children[c], sep="", file=file, append=TRUE)
  	}
    cat("\n", file=file, append=TRUE)
  } else {
    cat("incomplete\n", file=file, append=TRUE)
  }
}



cat.nodelist <- function(file=NULL) {
  cat("Analysis Structure as list of nodes\n", file=file)
  for (n in 1:length(TheAnalysisStructureNodes)) {
    cat.node(n, file)
  }
}





get.all.kind.and.ids <- function(A, kinds, ids) {
  if (is.list(A)) {
    if (is.character(A$kind)  &&  is.character(A$id)) {
      wk <- which(kinds == A$kind)
      wi <- integer(0)
      if (length(wk) > 0) {
        wi <- which(ids[wk] == A$id)
      }
      if (length(wi) == 0) {
        kinds <- c(kinds, A$kind)
        ids <- c(ids, A$id)
      }
    }
    for (i in 1:length(A)) {   
      ki <- get.all.kind.and.ids(A[[i]], kinds, ids)
      kinds <- ki$kinds
      ids <- ki$ids      
    }
  } 
  list(kinds=kinds, ids=ids)
}



########################################################################
######################internal functions ###############################
########################################################################


get.standard.nodes <- function(kind) {
	nodes <- NULL
	nodeIDs <- character(0)
	for (a in 1:get.nof.alignments()) {
		node <- TheAnalysisStructureNodes[[ ThePartitionIndices[a, kind] ]]
		if (!(node$id %in% nodeIDs)) {
			nodeIDs <- c(nodeIDs, node$id)
			nodes <- c(nodes, list(node))	
		}
	}
	nodes
}



make.new.node <- function(kind, id) {
  list(kind=kind, id=id, link=NULL, value=NULL, children=NULL)
}

make.link <- function(kind, id, link) {
  stopifnot(!is.null(link))
  list(kind=kind, id=id, link=link, value=NULL, children=NULL)
}

make.leaf <- function(kind, id, value) {
  stopifnot(is.atomic(value))
  stopifnot(!is.null(value))
  list(kind=kind, id=id, link=NULL, value=value, children=NULL)
}


add.new.node <- function(kind, id) {
  index <- look.for.kind.id(kind, id)
  if (index <= 0) {
    TheAnalysisStructureNodes <<- c(TheAnalysisStructureNodes, list(make.new.node(kind,  id)))
  } else {
    TheAnalysisStructureNodes <<- c(TheAnalysisStructureNodes, list(make.link(kind,  id, index)))
  }
  n <- length(TheAnalysisStructureNodes)
  if (index <= 0) {
    TheNodeHashTable[[paste0(kind, "__", id)]] <<- n
  }
  n
}


add.leaf <- function(kind, id, value) {
  index <- look.for.kind.id(kind, id)
  if (index <= 0) {
    TheAnalysisStructureNodes <<- c(TheAnalysisStructureNodes, list(make.leaf(kind,  id, value)))
  } else {
    TheAnalysisStructureNodes <<- c(TheAnalysisStructureNodes, list(make.link(kind,  id, index)))
  }
  n <- length(TheAnalysisStructureNodes)
  if (index <= 0) {
    TheNodeHashTable[[paste0(kind, "__", id)]] <<- n
  }
  n
}


fill.in.value <- function(kind, id, value) {
  idx <- look.for.kind.id(kind, id)
  stopifnot(idx > 0)
  TheAnalysisStructureNodes[[idx]] <<- 
    list(kind=TheAnalysisStructureNodes[[idx]]$kind, 
         id=TheAnalysisStructureNodes[[idx]]$id, 
         link=NULL, value=value, children=NULL) 
}


attach.children <- function(idx, children) {
  TheAnalysisStructureNodes[[idx]] <<- 
    list(kind=TheAnalysisStructureNodes[[idx]]$kind, 
         id=TheAnalysisStructureNodes[[idx]]$id, 
         link=NULL, value=NULL, children=children) 
}


index.from.kind.id <- function(kind, id) {
  n <- TheNodeHashTable[[paste0(kind, "__", id)]]
  stopifnot(!is.null(n))
  n
  
#   indices <- integer(0)
#   for (i in 1:length(TheAnalysisStructureNodes)) {
#     if (is.null(TheAnalysisStructureNodes[[i]]$link) &&
#           identical(TheAnalysisStructureNodes[[i]]$kind, kind)  &&
#           identical(TheAnalysisStructureNodes[[i]]$id, id)) {
#       indices <- c(indices, i)
#     }
#   }
#   if (length(indices)!=1) {
#     browser()
#   }
#   stopifnot(length(indices)==1)
#   n <- TheNodeHashTable[[paste0(kind, "__", id)]]
#   stopifnot(is.null(TheAnalysisStructureNodes[[n]]$link))
#   stopifnot(n == indices[1])
#   indices[1]
}


look.for.kind.id <- function(kind, id) {
  n <- TheNodeHashTable[[paste0(kind, "__", id)]]
  if (is.null(n)) { n <- -1 }
  n
  
  
#   indices <- integer(0)
#   for (i in 1:length(TheAnalysisStructureNodes)) {
#     if (is.null(TheAnalysisStructureNodes[[i]]$link) &&
#           identical(TheAnalysisStructureNodes[[i]]$kind, kind)  &&
#           identical(TheAnalysisStructureNodes[[i]]$id, id)) {
#       indices <- c(indices, i)
#     }
#   }
#   stopifnot(length(indices)<=1)
#   index <- -1
#   if (length(indices) == 1) {
#     index <- indices[1]
#   }
#   n <- TheNodeHashTable[[paste0(kind, "__", id)]]
#   if (index < 0) {
#     stopifnot(is.null(n))
#   } else {
#     stopifnot(n == index)
#   }
#   index
}



# Distinguishes between placeholder x and x with children.
bbas.real.kids <- function(x) {
  stopifnot(is.list(x))
  rk <- FALSE
  rk2 <- TRUE
  for (i in 3:length(x)) {
    if (is.null(x[[i]])) {
      rk2 <- FALSE
    } else {
      rk <- TRUE
    }
  }
  if (rk != rk2) {
  	browser()
    stop(paste0(bbas.kind.id(x), " is incomplete."))
  }
  rk
}





