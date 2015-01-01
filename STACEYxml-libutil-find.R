

nof.alignments <- function() {
  length(TheAnalysisStructure$alignment.table)
}



get.taxa.table <- function() {
  TheAnalysisStructure$taxa.table$taxa.table
}


get.minclusters <- function() {
  tt <- get.taxa.table()
  mcs <- tt[,2]
  unique(mcs)
}


get.bdc.model <- function() {
  # TODO not good that access first gtree for all
  TheAnalysisStructure$alignment.table$alignments[[1]]$gtree$prior$branching
}


get.bdcm.growthrate.fixed <- function() {
  (get.bdc.model()$growthrate$kind == "FixedValue")
}

get.bdcm.reldeath.fixed <- function() {
  (get.bdc.model()$reldeath$kind == "FixedValue")
}


get.bdcm.w.fixed <- function() {
  (get.bdc.model()$w$kind == "FixedValue")
}



get.smc.coalescent <- function() {
  # TODO not good that access first gtree for all
  TheAnalysisStructure$alignment.table$alignments[[1]]$gtree$prior$smccoal
}

get.smct.popsf.fixed <- function() {
  (get.smc.coalescent()$popSF$kind == "FixedValue")
}


AlignmentIndexToGtreeIndex <- NULL
TheGtreeModels <- NULL 

AlignmentIndexToClockIndex <- NULL
TheClockModels <- NULL

AlignmentIndexToSubstIndex <- NULL
TheSubstModels <- NULL 

AlignmentIndexToSiteHetIndex <- NULL
TheSiteHetModels <- NULL 




get.gtrees <- function() {
  if (is.null(TheGtreeModels)) {
    gtrees <- NULL
    gtreeids <- character(0)
    for (a in 1:nof.alignments()) {
      gte <- TheAnalysisStructure$alignment.table$alignments[[a]]$gtree
      if (length(which(gtreeids==gte$id)) == 0) {
        gtreeids <- c(gtreeids, gte$id)
        gtrees <- c(gtrees, list(gte))
      }
    }
    TheGtreeModels <<- gtrees
  }
  TheGtreeModels
}
  

g.of.alignment <- function(a) {
  if (is.null(AlignmentIndexToGtreeIndex)) {
    AlignmentIndexToGtreeIndex <- rep(-1, nof.alignments())
    gtrees <- get.gtrees()
    for (a in 1:nof.alignments()) {
      gte <- TheAnalysisStructure$alignment.table$alignments[[a]]$gtree
      for (g in 1:length(gtrees)) {
        if (gte$id == gtrees[[g]]$id) {
          AlignmentIndexToGtreeIndex[a] <- g
        }
      }
    }
  }
  AlignmentIndexToGtreeIndex[a]
}




get.clocks <- function() {
  if (is.null(TheClockModels)) {
    clocks <- NULL
    clockids <- character(0)
    for (a in 1:nof.alignments()) {
      clk <- TheAnalysisStructure$alignment.table$alignments[[a]]$clock
      if (length(which(clockids==clk$id)) == 0) {
        clockids <- c(clockids, clk$id)
        clocks <- c(clocks, list(clk))
      }
    }
    TheClockModels <<- clocks
  }
  TheClockModels
}

# TODO need c.of.gtree()
c.of.alignment <- function(a) {
  if (is.null(AlignmentIndexToClockIndex)) {
    AlignmentIndexToClockIndex <- rep(-1, nof.alignments())
    clocks <- get.clocks()
    for (a in 1:nof.alignments()) {
      clk <- TheAnalysisStructure$alignment.table$alignments[[a]]$clock
      for (g in 1:length(clocks)) {
        if (clk$id == clocks[[g]]$id) {
          AlignmentIndexToClockIndex[a] <- g
        }
      }
    }
  }
  AlignmentIndexToCSlockIndex[a]
}




get.substs <- function() {
  if (is.null(TheSubstModels)) {
    substs <- NULL
    substids <- character(0)
    for (a in 1:nof.alignments()) {
      sbt <- TheAnalysisStructure$alignment.table$alignments[[a]]$subst
      if (length(which(substids==sbt$id)) == 0) {
        substids <- c(substids, sbt$id)
        substs <- c(substs, list(sbt))
      }
    }
    TheSubstModels <<- substs
  }
  TheSubstModels
}


# TODO need b.of.gtree()
b.of.alignment <- function(a) {
  if (is.null(AlignmentIndexToSubstIndex)) {
    AlignmentIndexToSubstIndex <- rep(-1, nof.alignments())
    substs <- get.substs()
    for (a in 1:nof.alignments()) {
      sbt <- TheAnalysisStructure$alignment.table$alignments[[a]]$subst
      for (g in 1:length(substs)) {
        if (sbt$id == substs[[g]]$id) {
          AlignmentIndexToSubstIndex[a] <- g
        }
      }
    }
  }
  AlignmentIndexToSubstIndex[a]
}





get.sitehets <- function() {
  if (is.null(TheSiteHetModels)) {
    sitehets <- NULL
    sitehetids <- character(0)
    for (a in 1:nof.alignments()) {
      gte <- TheAnalysisStructure$alignment.table$alignments[[a]]$sitehet
      if (length(which(sitehetids==gte$id)) == 0) {
        sitehetids <- c(sitehetids, gte$id)
        sitehets <- c(sitehets, list(gte))
      }
    }
    TheSiteHetModels <<- sitehets
  }
  TheSiteHetModels
}


# TODO need h.of.gtree()

h.of.alignment <- function(a) {
  if (is.null(AlignmentIndexToSiteHetIndex)) {
    AlignmentIndexToSiteHetIndex <- rep(-1, nof.alignments())
    sitehets <- get.sitehets()
    for (a in 1:nof.alignments()) {
      gte <- TheAnalysisStructure$alignment.table$alignments[[a]]$sitehet
      for (g in 1:length(sitehets)) {
        if (gte$id == sitehets[[g]]$id) {
          AlignmentIndexToSiteHetIndex[a] <- g
        }
      }
    }
  }
  AlignmentIndexToSiteHetIndex[a]
}

